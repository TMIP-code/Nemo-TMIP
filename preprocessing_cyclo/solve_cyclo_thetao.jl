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
import Pardiso
using NonlinearSolve

#Adapted from TMIP-ACCESS solver scripts

########################################################################
# # # # # # # # # # # # # SOLVE THETAO # # # # # # # # # # # # # # # #
########################################################################


inputdir = ARGS[1] #path to input files
exp_tag = ARGS[2] #experiment identifier so that outputs do not overwrite
nprocs = 48

steps = collect(1:12)
Nsteps = length(steps)
δt = ustrip(s, 1yr / Nsteps)

########################################################################

volcello_ds = open_dataset(joinpath(inputdir, "volcello_.nc"))
areacello_ds = open_dataset(joinpath(inputdir, "areacello_.nc"))

# Build source variable for thetao
thetao_ds = open_dataset(joinpath(inputdir, "thetao_.nc"))
@time "building Temperatures" thetaos = [
    begin
            thetao = readcubedata(thetao_ds.thetao[month = At(month)])
	    thetao
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
thetao_mean = mean(thetaos)
src = Ω * thetao_mean[wet3D]

@time "initial state solve" u0 = solve(LinearProblem(M̄, src), solver, rtol = 1.0e-8).u
@show norm(M̄ * u0 - src) / norm(src)

temperatureinit3D = DimensionalData.rebuild(
    volcello_ds["volcello"];
    data = OceanTransportMatrixBuilder.as3D(u0, wet3D),
    dims = dims(volcello),
    metadata = Dict(
        "description" => "steady-state potential temperature stationary guess",
        "units" => "deg C",
    )
)

########################################################################
# Build seasonal source -- allow surface temperatures to vary as well as transport
@time "building seasonal source" src_cyclo = [
    begin
            thetao_ = thetaos[m]
	    src_ = Ω * thetao_[wet3D]
	    src_
        end
        for m in steps
]

#define stepping functions
function initstepprob(A, src)
    prob = LinearProblem(A, δt * src)
    return init(prob, solver, rtol = 1.0e-8)
end

function mystep!(du, u, p, m)
    prob = p.stepprob[m]
    prob.b = u .+ p.δt * p.src_cyclo[m] # xₘ₊₁ = Aₘ₊₁⁻¹ (xₘ + δt Ω x_srf)
    du .= solve!(prob).u
    return du
end
function jvpstep!(dv, v, p, m)
    prob = p.stepprob[m]
    prob.b = v # xₘ₊₁ = Aₘ₊₁⁻¹ (xₘ + δt 1)
    dv .= solve!(prob).u
    return dv
end
function steponeyear!(du, u, p)
    du .= u
    for m in eachindex(p.stepprob)
        mystep!(du, du, p, m)
    end
    return du
end
function jvponeyear!(dv, v, p)
    dv .= v
    for m in eachindex(p.stepprob)
        jvpstep!(dv, dv, p, m)
    end
    return dv
end
function G!(du, u, p)
    steponeyear!(du, u, p)
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
    stepprob::Vector{Any}
    src_cyclo::Vector{Vector{Float64}}
end

Base.copy(p::Params) = Params(p.δt, copy(stepprob), p.src_cyclo)

stepprob = [initstepprob(I + δt * Ms[m], src_cyclo[m]) for m in steps]
p = Params(δt, stepprob, src_cyclo)

nonlinearprob! = NonlinearProblem(f!, u0, p)


@info "solve cyclo-stationary state"
@time sol! = solve(nonlinearprob!, NewtonRaphson(linsolve = KrylovJL_GMRES(precs = precs, rtol = 1.0e-10)); show_trace = Val(true), reltol = Inf, abstol = 1.0e-10norm(u0, Inf));

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
                (m > 1) && mystep!(du, du, p, m)
                temperature3D = OceanTransportMatrixBuilder.as3D(du, wet3D)
                temperature4D = reshape(temperature3D, (size(wet3D)..., 1))
                axlist = (dims(volcello_ds["volcello"])..., dims(DimArray(ones(Nsteps), Ti(steps)))[1][m:m])
                temperature_YAXArray = rebuild(
                    volcello_ds["volcello"];
                    data = temperature4D,
                    dims = axlist,
                    metadata = Dict(
                        "origin" => "cyclo-stationary potential temperature by month",
                    )
                )
            end
            for m in eachindex(steps)
    )
)

temperaturemean3D = mean(cube4D, dims = Ti)

arrays = Dict(:thetao => temperaturemean3D, :thetao_init => temperatureinit3D, :lat => volcello_ds.lat, :lon => volcello_ds.lon)
ds = Dataset(; volcello_ds.properties, arrays...)

# Save temperature3D to netCDF file
outputfile = joinpath(inputdir, "mean_cyclomon_thetao$exp_tag.nc")
@info "Saving potential temperature as netCDF file:\n  $(outputfile)"
savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)

