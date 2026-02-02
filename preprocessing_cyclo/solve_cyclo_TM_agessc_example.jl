using OceanTransportMatrixBuilder
using NetCDF
using YAXArrays
using DataFrames
using DimensionalData
using SparseArrays
using LinearAlgebra
using Unitful
using Unitful: s, yr, d
using Statistics
using Format
using Dates
using FileIO
using LinearSolve
using IncompleteLU
import Pardiso
using NonlinearSolve
using MPI
import MUMPS

#Adapted from TMIP-ACCESS solver scripts

########################################################################
# # # # # # # # # # # # # SOLVE AGESSC # # # # # # # # # # # # # # # #
########################################################################


inputdir = ARGS[1] #path to input files
exp_tag = ARGS[2] #experiment identifier so that outputs do not overwrite
nprocs = 24

steps = collect(1:12)
Nsteps = length(steps)
δt = ustrip(s, 1yr / Nsteps)

########################################################################

volcello_ds = open_dataset(joinpath(inputdir, "volcello_.nc"))
areacello_ds = open_dataset(joinpath(inputdir, "areacello_.nc"))

# Build source variable for agessc
agessc_ds = open_dataset(joinpath(inputdir, "agessc_.nc"))
@time "building Ages" agesscs = [
    begin
            agessc = readcubedata(agessc_ds.agessc[month = At(month)])
	    agessc
        end
        for month in steps
]

# Load fixed variables in memory
areacello = readcubedata(areacello_ds.areacello)
volcello = readcubedata(volcello_ds.volcello)
lon = readcubedata(volcello_ds.lon)
lat = readcubedata(volcello_ds.lat)
lev = volcello_ds.lev
lon_vertices = volcello_ds.lon_verticies
lat_vertices = volcello_ds.lat_verticies

gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
(; lon_vertices, lat_vertices, v3D) = gridmetrics

indices = makeindices(v3D)
(; wet3D, L, Lwet, N, Lwet3D, C) = indices

(; zt) = gridmetrics
########################################################################

#set surface source region masks

issrf = let #mask of 'is surface'
    issrf3D = falses(size(wet3D))
    issrf3D[:, :, 1] .= true
    issrf3D[wet3D]
end

Ω = sparse(Diagonal(Float64.(issrf)))

########################################################################

# The equation for conservative θ or S is
#   ∂ₜx + T x = Ω (x_srf - x)
# Applying Backward Euler time step gives
#   (I + Δt M) xₖ₊₁ = xₖ + Δt Ω x_srf
# where M = T + Ω



# Build matrices -- grab matrixes with exp_tag
@time "building Ms" Ms = [
    begin
            inputfile = joinpath(inputdir, "cyclo_matrix_month$m$exp_tag.jld2")
            @info "Loading matrices + metrics as $inputfile"
            T = load(inputfile)["T"]
            T + Ω
        end
        for m in steps
]


#######################################################################

# Preconditioner -- precondition the cyclostationary matrix on the stationary case
M̄ = mean(Ms) #
Δt = sum(δt for _ in steps)

#This is the preconditioner from Benoit's TMIP-ACCESS scripts
struct CycloPreconditioner
    prob
end
Base.eltype(::CycloPreconditioner) = Float64
function LinearAlgebra.ldiv!(Pl::CycloPreconditioner, x::AbstractVector)
    @info "applying Pl"
    Pl.prob.b = x
    solve!(Pl.prob)
    x .= Pl.prob.u .- x # Note the -x (following Bardin et al)
    return x
end
function LinearAlgebra.ldiv!(y::AbstractVector, Pl::CycloPreconditioner, x::AbstractVector)
    Pl.prob.b = x
    solve!(Pl.prob)
    y .= Pl.prob.u .- x # Note the -x (following Bardin et al)
    return y
end

matrix_type = Pardiso.REAL_SYM
solver = MKLPardisoIterate(; nprocs, matrix_type)

Plprob = LinearProblem(-Δt * M̄, ones(N))  # following Bardin et al. (M -> -M though)
Plprob = init(Plprob, solver, rtol = 1.0e-8)
Pl = CycloPreconditioner(Plprob)
Pr = I
precs = Returns((Pl, Pr))

########################################################################

# use static (annual mean) solution as initial guess

#The MUMPS solver handles the centered scheme case well
MPI.Init()
@show norm(M̄ * u0 - ones(N)) / norm(ones(N))
@show norm(M̄ * u0 - src) / norm(src)
MPI.Finalize()
u0 = vec(u0)

ageinit3D = DimensionalData.rebuild(
    volcello_ds["volcello"];
    data = ustrip.(yr, OceanTransportMatrixBuilder.as3D(u0, wet3D) * s),
    dims = dims(volcello),
    metadata = Dict(
        "description" => "steady-state ideal mean age stationary guess",
        "units" => "yr",
    )
)

########################################################################


#define stepping functions
function initstepprob(A)
    prob = LinearProblem(A, ones(N))
    return init(prob, solver, rtol = 1.0e-8)
end

stepprob = [initstepprob(I + δt * M) for M in Ms]

function stepforwardonemonth!(du, u, p, m)
    prob = stepprob[m]
    prob.b = u .+ p.δt #prob.b = u .+ p.δt # xₘ₊₁ = Aₘ₊₁⁻¹ (xₘ + δt 1) # CHECK m index is not off by 1
    du .= solve!(prob).u
    return du
end
function jvpstep!(dv, v, p, m)
    prob = stepprob[m]
    prob.b = v # xₘ₊₁ = Aₘ₊₁⁻¹ (xₘ + δt 1) # CHECK m index is not off by 1
    dv .= solve!(prob).u
    return dv
end
function stepforwardoneyear!(du, u, p)
    du .= u
    for m in steps
        stepforwardonemonth!(du, du, p, m)
    end
    return du
end
function jvponeyear!(dv, v, p)
    dv .= v
    for m in steps
        jvpstep!(dv, dv, p, m)
    end
    return dv
end
function G!(du, u, p)
    stepforwardoneyear!(du, u, p)
    du .-= u
    return du
end
function jvp!(dv, v, u, p)
    jvponeyear!(dv, v, p)
    dv .-= v
    return dv
end

########################################################################
#set up nonlinear problem

f! = NonlinearFunction(G!; jvp = jvp!)

struct Params
    δt::Float64
end

Base.copy(p::Params) = Params(p.δt)

p = Params(δt)

nonlinearprob! = NonlinearProblem(f!, u0, p)

@info "solve cyclo-stationary state"
@time sol! = solve(nonlinearprob!, NewtonRaphson(linsolve = KrylovJL_GMRES(precs = precs, rtol = 1.0e-10)); show_trace = Val(true), reltol = Inf, abstol = 1.0e-8norm(u0, Inf));

@info "Check the RMS drift, should be order 10⁻¹¹‰ (1e-11 per thousands)"
du = deepcopy(u0)
@show norm(G!(du, sol!.u, p), Inf) / norm(sol!.u, Inf) |> u"permille"

########################################################################
#Save solution

du = sol!.u
cube4D = reduce(
    (a, b) -> cat(a, b, dims = Ti),
    (
        begin
                (m > 1) && stepforwardonemonth!(du, du, p, m)
		duinyr = ustrip.(yr, du .* s)
                age3D = OceanTransportMatrixBuilder.as3D(duinyr, wet3D)
                age4D = reshape(age3D, (size(wet3D)..., 1))
                axlist = (dims(volcello_ds["volcello"])..., dims(DimArray(ones(Nsteps), Ti(steps)))[1][m:m])
                age_YAXArray = rebuild(
                    volcello_ds["volcello"];
                    data = age4D,
                    dims = axlist,
                    metadata = Dict(
                        "origin" => "cyclo-stationary ideal mean age by month",
                    )
                )
            end
            for m in eachindex(steps)
    )
)

agemean3D = mean(cube4D, dims = Ti)

arrays = Dict(:agessc => agemean3D, :agessc_init => ageinit3D, :lat => volcello_ds.lat, :lon => volcello_ds.lon)
ds = Dataset(; volcello_ds.properties, arrays...)

# Save age3D to netCDF file
outputfile = joinpath(inputdir, "mean_cyclomon_agessc$exp_tag.nc")
@info "Saving ideal age as netCDF file:\n  $(outputfile)"
savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)

