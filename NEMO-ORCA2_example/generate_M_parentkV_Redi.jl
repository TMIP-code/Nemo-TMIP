using OceanTransportMatrixBuilder
using NetCDF
using YAXArrays
using DataFrames
using DimensionalData
using SparseArrays
using LinearAlgebra
using Unitful
using Unitful: s, yr
using Statistics
using NaNStatistics
using Format
using Dates
using FileIO
include("../OceanTransportMatrixBuilder_extrafunctions/matrixbuilding_extras_Redi.jl")
include("../OceanTransportMatrixBuilder_extrafunctions/matrixbuilding_extras_parentkV.jl")


############### Build the transport matrix #############################
data_path = "../data/processed/NEMO-ORCA2/"
months = collect(1:12)

umo_ds = open_dataset(joinpath(data_path, "umo.nc"))
vmo_ds = open_dataset(joinpath(data_path, "vmo.nc"))
mlotst_ds = open_dataset(joinpath(data_path, "maxmld.nc"))
rho_ds = open_dataset(joinpath(data_path, "rho.nc"))
volcello_ds = open_dataset(joinpath(data_path, "volcello.nc"))
areacello_ds = open_dataset(joinpath(data_path, "areacello.nc"))
avt_ds = open_dataset(joinpath(data_path, "kV.nc"))


#load fixed variables into memory
areacello = readcubedata(areacello_ds.areacello)
volcello = readcubedata(volcello_ds.volcello)
lon = readcubedata(volcello_ds.lon)
lat = readcubedata(volcello_ds.lat)
lev = volcello_ds.lev # depth of cell centers
lon_vertices = areacello_ds.lon_verticies # cell vertices
lat_vertices = areacello_ds.lat_verticies # cell vertices


gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
(; lon_vertices, lat_vertices) = gridmetrics
indices = makeindices(gridmetrics.v3D)


# parameter values
ρ = 1035.0       # density (kg/m^3)
κV = 1.4e-7  #provide a background field for stability -- molecular floor
κH = 375 #nominal isopycnal diffusivity (m^2/s)
Ld =  200.e+3 #nominal grid spacing (m) for scaling isopycnal diffusivity
κVML = 1. # mixed layer diffusivity (m^2/s) NOT USED
upwind=false #use centered scheme

#calculate "diffusive velocity scale" to feed to redi_build_operators
Ud = 2 * κH / Ld



for month in months

    # Load variables in memory
    mlotst = readcubedata(mlotst_ds.mlotstmax[month = At(month)])
    umo = readcubedata(umo_ds.umo[month = At(month)])
    vmo = readcubedata(vmo_ds.vmo[month = At(month)])
    rho3D = readcubedata(rho_ds.rho2[month = At(month)])
    avt = readcubedata(avt_ds.avt[month = At(month)])
    
    # Also remove missings in umo, vmo, rho3D, avt
    umo = replace(umo, missing => 0)
    vmo = replace(vmo, missing => 0)
    rho3D = replace(rho3D, missing => 0)
    avt = replace(avt, missing => 0)
    avt = replace(avt, NaN => 0.)

    ϕ = facefluxesfrommasstransport(; umo, vmo, gridmetrics, indices)


    # Make the transport matrix
    (; Tadv, TκH, TκVML, TκVdeep)  = transportmatrix(; ϕ, mlotst, gridmetrics, indices, ρ, κH, κVML, κVdeep=κV, upwind)

    (; TκHredi, TκVredi, TκVcross) = redi_build_operators(; gridmetrics, indices, rho3D, mlotst, ldf_c2d = true, rn_Ud = Ud)

    Tκavt = buildTκVavt(; gridmetrics, indices, avt)

    # add up advection (Tadv), parent model kV (Tκavt), background molecular kV (TκVdeep)
    # and all components of the Redi diffusion (TκHredi + TκVredi + TκVcross)
    T = Tadv + Tκavt + TκVdeep + TκHredi + TκVredi + TκVcross
    
    outputfile = joinpath(data_path, "M_$month.jld2")
    @info "Saving matrix as $outputfile"
    save(
		outputfile,
        	Dict(
            	"T" => T,
        	)
    )

end