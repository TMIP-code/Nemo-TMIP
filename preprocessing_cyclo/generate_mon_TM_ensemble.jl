using OceanTransportMatrixBuilder
using NetCDF
using YAXArrays
using DataFrames
using DimensionalData
using SparseArrays
using LinearAlgebra
using Unitful
using Unitful: s, yr
using NaNStatistics
using Format
using Dates
using FileIO

data_path = ARGS[1]
exp_tag = ARGS[2]
κH = parse(Float64, ARGS[3])
κVdeep = parse(Float64, ARGS[4])
############### Build the transport matrix #############################
months = collect(1:12)

umo_ds = open_dataset(joinpath(data_path, "umo_.nc"))
vmo_ds = open_dataset(joinpath(data_path, "vmo_.nc"))
mlotst_ds = open_dataset(joinpath(data_path, "mlotst_.nc"))
volcello_ds = open_dataset(joinpath(data_path, "volcello_.nc"))
areacello_ds = open_dataset(joinpath(data_path, "areacello_.nc"))

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

# Some parameter values  -- use the parameters from Pasquier et al. 2025
ρ = 1035.0     # density (kg/m^3)
κVML = 100.    # mixed-layer vertical diffusivity (m^2/s)
upwind = false

for month in months

    # Load variables in memory
    mlotst = readcubedata(mlotst_ds.mlotst[month = At(month)])
    umo = readcubedata(umo_ds.umo[month = At(month)])
    vmo = readcubedata(vmo_ds.vmo[month = At(month)])

    # Also remove missings in umo and vmo
    umo = replace(umo, missing => 0)
    vmo = replace(vmo, missing => 0)

    ϕ = facefluxesfrommasstransport(; umo, vmo, gridmetrics, indices)

    # Make the transport matrix (should take a few seconds)
    (; T) = transportmatrix(; ϕ, mlotst, gridmetrics, indices, ρ, κH, κVML, κVdeep, upwind)

    
    outputfile = joinpath(data_path, "cyclo_matrix_month$month$exp_tag.jld2")
    @info "Saving matrix as $outputfile"
    save(
        outputfile,
        Dict(
            "T" => T,
        )
    )

end