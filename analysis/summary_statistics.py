import xarray as xr
import numpy as np

#this script compares the correlation and biases of T,S for different matrix experiments
# as well as moments of the global age distribution

path = '../data/processed/cyclo_stationary/IPSL-CM6A-LR_r1i2p1f1_monthly/'
exp_tags = ["cntrl"] #["_cntrl", "_exp1", "_exp2", "_exp3", "_exp4", "_exp5", "_exp6", "_exp7", "_exp8"]

#online fields
time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
DS1 = xr.open_dataset(path + 'thetao_.nc', decode_times=time_coder)
DS1['so'] =  xr.open_dataset(path + 'so_.nc', decode_times=time_coder).so
#uncomment if parent model provides agessc output
#DS1['agessc'] =  xr.open_dataset(path + 'agessc_.nc', decode_times=time_coder).agessc
DS1 = DS1.mean('month')


for exp in exp_tags:
    print(f'########################## Stats for experiment tag: {exp} #########################')

    #matrix fields
    DS2 = xr.open_dataset(path + f'mean_cyclomon_thetao{exp}.nc')
    DS2['so'] = xr.open_dataset(path + f'mean_cyclomon_so{exp}.nc').so
    DS2['so_init'] = xr.open_dataset(path + f'mean_cyclomon_so{exp}.nc').so_init
    DS2['agessc'] = xr.open_dataset(path + f'mean_cyclomon_agessc{exp}.nc').agessc
    DS2['agessc_init'] = xr.open_dataset(path + f'mean_cyclomon_agessc{exp}.nc').agessc_init    
    DS2 = DS2.mean('Ti').set_coords(['lat', 'lon'])

    vol = xr.open_dataset(path+'volcello_.nc', decode_times=time_coder).volcello
    area = xr.open_dataset(path+'areacello_.nc', decode_times=time_coder).areacello
    thk = vol/area

    #Upper Ocean (100-500m)
    depth_slice = slice(100, 500)

    thetao_corr = xr.corr(DS2.thetao.sel(lev = depth_slice), DS1.thetao.sel(lev = depth_slice), weights = vol.fillna(0.))
    thetao_corr0 = xr.corr(DS2.thetao_init.sel(lev = depth_slice), DS1.thetao.sel(lev = depth_slice), weights = vol.fillna(0.))
    print(f"T correlation (100-500m): {thetao_corr.values:.3f} (initial guess {thetao_corr0.values:.3f})")

    so_corr = xr.corr(DS2.so.sel(lev = depth_slice), DS1.so.sel(lev = depth_slice), weights = vol.fillna(0.))
    so_corr0 = xr.corr(DS2.so_init.sel(lev = depth_slice), DS1.so.sel(lev = depth_slice), weights = vol.fillna(0.))
    print(f"S correlation (100-500m): {so_corr.values:.3f} (initial guess {so_corr0.values:.3f})")

    thetao_bias = (DS2.thetao.sel(lev = depth_slice) - DS1.thetao.sel(lev = depth_slice)).weighted(vol.fillna(0.)).mean()
    thetao_bias0 = (DS2.thetao_init.sel(lev = depth_slice) - DS1.thetao.sel(lev = depth_slice)).weighted(vol.fillna(0.)).mean()
    print(f"T bias (100-500m): {thetao_bias.values:.3f} (initial guess {thetao_bias0.values:.3f})")

    so_bias = (DS2.so.sel(lev = depth_slice) - DS1.so.sel(lev = depth_slice)).weighted(vol.fillna(0.)).mean()
    so_bias0 = (DS2.so_init.sel(lev = depth_slice) - DS1.so.sel(lev = depth_slice)).weighted(vol.fillna(0.)).mean()
    print(f"S bias (100-500m): {so_bias.values:.3f} (initial guess {so_bias0.values:.3f})")

    #Full depth
    thetao_corr = xr.corr(DS2.thetao, DS1.thetao, weights = vol.fillna(0.))
    thetao_corr0 = xr.corr(DS2.thetao_init, DS1.thetao, weights = vol.fillna(0.))
    print(f"T correlation (full depth): {thetao_corr.values:.3f} (initial guess {thetao_corr0.values:.3f})")

    so_corr = xr.corr(DS2.so, DS1.so, weights = vol.fillna(0.))
    so_corr0 = xr.corr(DS2.so_init, DS1.so, weights = vol.fillna(0.))
    print(f"S correlation (full depth): {so_corr.values:.3f} (initial guess {so_corr0.values:.3f})")

    thetao_bias = (DS2.thetao - DS1.thetao).weighted(vol.fillna(0.)).mean()
    thetao_bias0 = (DS2.thetao_init - DS1.thetao).weighted(vol.fillna(0.)).mean()
    print(f"T bias (full depth): {thetao_bias.values:.3f} (initial guess {thetao_bias0.values:.3f})")

    
    so_bias = (DS2.so - DS1.so).weighted(vol.fillna(0.)).mean()
    so_bias0 = (DS2.so_init - DS1.so).weighted(vol.fillna(0.)).mean()
    print(f"S bias (full depth): {so_bias.values:.3f} (initial guess {so_bias0.values:.3f})")

    mean_age = DS2.agessc.weighted(vol.fillna(0.)).mean()
    mean_age0 = DS2.agessc_init.weighted(vol.fillna(0.)).mean()
    print(f"Mean Age (full depth): {mean_age.values:.4f} (initial guess {mean_age0.values:.3f})")

    std_age = DS2.agessc.weighted(vol.fillna(0.)).std()
    std_age0 = DS2.agessc_init.weighted(vol.fillna(0.)).std()
    print(f"Std Age (full depth): {std_age.values:.4f} (initial guess {std_age0.values:.3f})")

    max_age = DS2.agessc.max()
    max_age0 = DS2.agessc_init.max()
    print(f"Max Age (full depth): {max_age.values:.4f} (initial guess {max_age0.values:.3f})")
