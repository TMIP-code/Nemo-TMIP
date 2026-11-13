import xarray as xr
import numpy as np

#process and standardize the data to be used as inputs for the OceanTransportMatrixBuilder

model = 'IPSL-CM6A-LR'
member = 'r1i2p1f1'
path = f'../data/processed/cyclo_stationary/{model}_{member}_monthly/'

#default is monthly climatology
#flag to output seasonal climatology (DJF, MAM, JJA, SON), instead
SSN = False

#compute density from temperature and salinity
compute_rho=False

#standardize the coordinate and dimension names
vrs_to_rename = {'olevel':'lev',
                 'nav_lat':'lat', 'nav_lon':'lon',
                 'bounds_nav_lon': 'lon_verticies', 'bounds_nav_lat':'lat_verticies',                 
                 'nav_lon_bnds': 'lon_verticies', 'nav_lat_bnds':'lat_verticies',
                 'nvertex':'vertices',
                 'olevel_bnds':'lev_bounds'}

time_vrs = ['time', 'time_bnds']

time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    
print('Clean datasets for {model} {member}')
vrs = ['vo', 'uo', 'areacello', 'mlotst', 'umo', 'vmo', 'thetao', 'so', 'agessc']
for vr in vrs:
    ds = xr.open_dataset(path + f'{vr}.nc', decode_times=time_coder)
    for v, v_ in vrs_to_rename.items():
        if v in list(ds.variables):
            ds = ds.rename({v:v_})
            ds = ds.assign_coords({v_:ds[v_]})
    ds[vr] = ds[vr].astype('float64')

    if 'time' in list(ds.variables):
        if vr in ['volcello']:
            ds = ds.mean('time')
        else:
            if SSN:
                ds = ds.groupby(ds.time.dt.season).mean('time')
                ds.mean('season').squeeze().to_netcdf(path + f'{vr}_mean.nc')
            else:
                ds = ds.groupby(ds.time.dt.month).mean('time')
        
        if 'time' in ds.encoding.get('unlimited_dims', []):
            ds.encoding['unlimited_dims'].remove('time')

    #remove overlapping longitude columns and relabel longitude 0-360
    ds = ds.isel(x =slice(1, 361))
    ds['lon'] = np.mod(ds.lon, 360)
    
    print(f'Saving {vr}...')
    ds.squeeze().to_netcdf(path + f'{vr}_.nc')


#construct volume from lev_bounds -- could be improved by incorporating info about ocean depth!!!
ds_example =  xr.open_dataset(path + 'thetao_.nc', decode_times=time_coder).mean('month')
ds_area = xr.open_dataset(path + 'areacello_.nc', decode_times=time_coder)
DZ = ds_example.lev_bounds.isel(bnds=1) - ds_example.lev_bounds.isel(bnds=0)
ds_example['volcello'] = (ds_area.areacello*DZ).where(~np.isnan(ds_example.thetao)).transpose('lev', 'y', 'x')
ds_example = ds_example.drop_vars(['thetao'])
print("Saving volcello...")
ds_example.to_netcdf(path + "volcello_.nc")

#compute density
if compute_rho:
    from fastjmd95 import rho
    SA = xr.open_dataset(path + 'so_.nc', decode_times=time_coder).so
    PT = xr.open_dataset(path + 'thetao_.nc', decode_times=time_coder).thetao
    rho2 = rho(SA, PT, 2000.)
    ds_rho = rho2.to_dataset(name = 'rho2')
    print('Saving rho_2...')
    ds_rho.to_netcdf(path + 'rho_.nc')
