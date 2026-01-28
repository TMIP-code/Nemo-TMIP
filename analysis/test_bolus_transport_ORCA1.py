import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy as cp

#Test whether umo/vmo include contributions from bolus velocities
## by comparing with vo/uo

model = 'IPSL-CM6A-LR'
member = 'r1i2p1f1'

path = f'../data/processed/cyclo_stationary/{model}_{member}_monthly/'

ds_uo = xr.open_dataset(path + 'uo_.nc').mean('month')
ds_vo = xr.open_dataset(path + 'vo_.nc').mean('month')

ds_umo = xr.open_dataset(path + 'umo_.nc').mean('month')
ds_vmo = xr.open_dataset(path + 'vmo_.nc').mean('month')

ds_vol = xr.open_dataset(path + 'volcello_.nc')
ds_area = xr.open_dataset(path + 'areacello_.nc')

#read in grid metrics for the ORCA1 tripolar grid
ds_orca = xr.open_dataset("mesh_mask_orca1.nc4").mean('t')

#generate supergrid and standardize names. I think the Nemo grid left handed?
u_grid_names = {'x':'x_l', 'y':'y_c', 'lon':'lon_u', 'lat': 'lat_u', 'lev':'z'}
v_grid_names = {'x':'x_c', 'y':'y_l', 'lon':'lon_v', 'lat': 'lat_v', 'lev':'z'}
t_grid_names = {'x':'x_c', 'y':'y_c', 'lon':'lon_c', 'lat': 'lat_c', 'lev':'z'}

ds_uo = ds_uo.rename(u_grid_names)[['uo']]
ds_umo = ds_umo.rename(u_grid_names)[['umo']]
ds_vo = ds_vo.rename(v_grid_names)[['vo']]
ds_vmo = ds_vmo.rename(v_grid_names)[['vmo']]
ds_vol = ds_vol.rename(t_grid_names)[['volcello']]
ds_area = ds_area.rename({'x':'x_c', 'y':'y_c', 'lon':'lon_c', 'lat': 'lat_c'})[['areacello']]

DS = xr.merge([ds_uo, ds_umo, ds_vo, ds_vmo, ds_vol, ds_area], compat = 'no_conflicts')

#extract grid dimensions
dx_t = ds_orca.e1t.rename({'x':'x_c', 'y':'y_c'})
dy_t = ds_orca.e2t.rename({'x':'x_c', 'y':'y_c'})
dz_t = ds_orca.e3t_0.rename({'x':'x_c', 'y':'y_c'})

dx_u = ds_orca.e1u.rename({'x':'x_l', 'y':'y_c'})
dy_u = ds_orca.e2u.rename({'x':'x_l', 'y':'y_c'})
dz_u = ds_orca.e3u_0.rename({'x':'x_l', 'y':'y_c'})

dx_v = ds_orca.e1v.rename({'x':'x_c', 'y':'y_l'})
dy_v = ds_orca.e2v.rename({'x':'x_c', 'y':'y_l'})
dz_v = ds_orca.e3v_0.rename({'x':'x_c', 'y':'y_l'})

DS['dx_t'] = dx_t.where(dx_t > 0.)
DS['dy_t'] = dy_t.where(dy_t > 0.)
DS['dz_t'] = dz_t.where(dx_t > 0.)

DS['dx_u'] = dx_u.where(dx_u > 0.)
DS['dy_u'] = dy_u.where(dy_u > 0.)
DS['dz_u'] = dz_u.where(dy_u > 0.)

DS['dx_v'] = dx_v.where(dx_v > 0.)
DS['dy_v'] = dy_v.where(dy_v > 0.)
DS['dz_v'] = dz_v.where(dy_v > 0.)

rho_0 = 1030. #the result isn't too sensitive here
mo_to_Sv = 1 / (1e6 * rho_0)
umo_from_uo = DS.uo * DS.dy_u * DS.dz_u * rho_0
vmo_from_vo = DS.vo * DS.dx_v * DS.dz_v * rho_0

#I don't have a reference for how close they would be if they included the same processes
#but I'm not sure this is a super useful metrics
rel_err_uo = np.abs((umo_from_uo - DS.umo)/DS.umo).weighted(DS.volcello.fillna(0.)).mean() * 100.
rel_err_vo = np.abs((vmo_from_vo - DS.vmo)/DS.vmo).weighted(DS.volcello.fillna(0.)).mean() * 100.
print(f'Weighted mean relative error for uo: {rel_err_uo.values:.3f} %')
print(f'Weighted mean relative error for vo: {rel_err_vo.values:.3f} %')


#try a more aggregated metric -- deep overturning at 40S
MOC = DS.vmo.sum('x_c').cumsum('z')*mo_to_Sv
MOC_from_vo = vmo_from_vo.sum('x_c').cumsum('z')*mo_to_Sv
#40S at y_l=126
MOC_ind = MOC.sel(y_l=126).sel(z = slice(2000, 5000)).min('z')
MOC_from_vo_ind = MOC_from_vo.sel(y_l=126).sel(z = slice(2000, 5000)).min('z')
print(f"deep overturning cell at 40S from vmo is: {MOC_ind.values:.3f} Sv")
print(f"deep overturning cell at 40S from vo is:  {MOC_from_vo_ind.values:.3f} Sv")

#Plot the bolus component
## if the bolus component is not included in umo/vmo
## then this should just look like noise

fig, axes = plt.subplots(2,3, figsize = [20, 10], subplot_kw={'projection':cp.crs.Robinson(205)})
cbar_kwargs = {'orientation':'horizontal', 'shrink':0.5, 'pad':0.05, 'label':'[Sv]', 'ticks':[-.05, 0., 0.05]}
kwargs = {'cmap':'RdBu_r',  'vmax':.05, 'vmin':-.05, 'levels':21, 'transform':cp.crs.PlateCarree(), 'cbar_kwargs':cbar_kwargs}

for ax in axes.flatten():
    ax.set_global()
    ax.add_feature(cp.feature.LAND, zorder=10)
    
for i, depths in enumerate([slice(10, 100), slice(300, 500), slice(3000, 5000)]):

    umo_from_uo_2d = (umo_from_uo.sel(z = depths).sum('z') * mo_to_Sv).rename({'x_l':'x', 'y_c':'y', 'lat_u':'lat', 'lon_u':'lon'})
    umo_2d = (DS.umo.sel(z = depths).sum('z') * mo_to_Sv).rename({'x_l':'x', 'y_c':'y', 'lat_u':'lat', 'lon_u':'lon'})
    diff_umo = umo_from_uo_2d - umo_2d

    vmo_from_vo_2d = (vmo_from_vo.sel(z = depths).sum('z') *mo_to_Sv).rename({'x_c':'x', 'y_l':'y', 'lat_v':'lat', 'lon_v':'lon'})
    vmo_2d = (DS.vmo.sel(z = depths).sum('z') * mo_to_Sv).rename({'x_c':'x', 'y_l':'y', 'lat_v':'lat', 'lon_v':'lon'})
    diff_vmo = vmo_from_vo_2d - vmo_2d

    ax = axes[0,i]
    diff_umo.plot(ax=ax, x='lon', y='lat', **kwargs)
    
    ax = axes[1,i]
    diff_vmo.plot(ax=ax, x='lon', y='lat', **kwargs)

axes[0,0].set_title('10-100 m [U]')
axes[1,0].set_title('10-100 m [V]')
axes[0,1].set_title('300-500 m [U]')
axes[1,1].set_title('300-500 m [V]')
axes[0,2].set_title('3000-5000 m [U]')
axes[1,2].set_title('3000-5000 m [V]')

plt.show()
