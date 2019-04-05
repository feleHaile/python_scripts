'''15/3/2019 Save JRA55 energy tendencies and MSE components in one file
'''
import xarray as xr
import sh
import numpy as np

def pentad_mean_climatology(data, years):  # Function to get pentad of year
    pentad_years = np.array([])
    for year in years:
        data_year = data.sel(time=str(year))
        if len(data_year.time)==366:
            pentad = np.repeat(np.arange(1., 74.), 5)
            pentad = np.insert(pentad, 10, 2)    
        else:
            pentad = np.repeat(np.arange(1., 74.), 5)    
        pentad_years = np.concatenate((pentad_years, pentad))
        
    data = data.assign_coords(pentad = ('time', pentad_years))
    
    data_pentads = data.groupby('pentad').mean(('time'))
    
    return data_pentads

    
ds = xr.open_dataset('/disca/share/rg419/JRA_55/heating_terms.nc')

#radiative heating
#data_swhr = xr.open_dataset('/disca/share/rg419/JRA_55/swhr_daily/atmos_daily_together_2_5.nc', chunks={'time': 30});     data_swhr = data_swhr['var250'].load().loc['1979-01':'2016-12'].sel(plev=85000.)
#data_swhr = pentad_mean_climatology(data_swhr, np.arange(1979,2017))
#data_lwhr = xr.open_dataset('/disca/share/rg419/JRA_55/lwhr_daily/atmos_daily_together_2_5.nc', chunks={'time': 30});     data_lwhr = data_lwhr['var251'].load().loc['1979-01':'2016-12'].sel(plev=85000.)
#data_lwhr = pentad_mean_climatology(data_lwhr, np.arange(1979,2017))

#ds = xr.Dataset({'swhr': (['pentad','lat', 'lon'], data_swhr),
#                 'lwhr':  (['pentad', 'lat', 'lon'], data_lwhr)},
#                     coords={'pentad': ('pentad', data_swhr.pentad),
#                               'lat': ('lat', data_swhr.lat),
#                               'lon': ('lon', data_swhr.lon)})

#print(ds)
    
# other heating
#data_cnvhr = xr.open_dataset('/disca/share/rg419/JRA_55/cnvhr_daily/atmos_daily_together_2_5.nc', chunks={'time': 30});     data_cnvhr = data_cnvhr['var242'].load().loc['1979-01':'2016-12'].sel(plev=85000.)
#ds['cnvhr'] = pentad_mean_climatology(data_cnvhr, np.arange(1979,2017))
#data_lrghr = xr.open_dataset('/disca/share/rg419/JRA_55/lrghr_daily/atmos_daily_together_2_5.nc', chunks={'time': 30});     data_lrghr = data_lrghr['var241'].load().loc['1979-01':'2016-12'].sel(plev=85000.)
#ds['lrghr'] = pentad_mean_climatology(data_lrghr, np.arange(1979,2017))    
#data_vdfhr = xr.open_dataset('/disca/share/rg419/JRA_55/vdfhr_daily/atmos_daily_together_2_5.nc', chunks={'time': 30});     data_vdfhr = data_vdfhr['var246'].load().loc['1979-01':'2016-12'].sel(plev=85000.)
#ds['vdfhr'] = pentad_mean_climatology(data_vdfhr, np.arange(1979,2017))
#data_adhr = xr.open_dataset('/disca/share/rg419/JRA_55/adhr_daily/atmos_daily_together_2_5.nc', chunks={'time': 30});     data_adhr = data_adhr['var222'].load().loc['1979-01':'2016-12'].sel(plev=85000.)
#ds['adhr'] = pentad_mean_climatology(data_adhr, np.arange(1979,2017))

#print(ds)

# moistening
#data_cnvmr = xr.open_dataset('/disca/share/rg419/JRA_55/cnvmr_daily/atmos_daily_together_2_5.nc', chunks={'time': 30});     data_cnvmr = data_cnvmr['var243'].load().loc['1979-01':'2016-12'].sel(plev=85000.)
#ds['cnvmr'] = pentad_mean_climatology(data_cnvmr, np.arange(1979,2017))
#data_lrgmr = xr.open_dataset('/disca/share/rg419/JRA_55/lrgmr_daily/atmos_daily_together_2_5.nc', chunks={'time': 30});     data_lrgmr = data_lrgmr['var253'].load().loc['1979-01':'2016-12'].sel(plev=85000.)
#ds['lrgmr'] = pentad_mean_climatology(data_lrgmr, np.arange(1979,2017))    
#data_vdfmr = xr.open_dataset('/disca/share/rg419/JRA_55/vdfmr_daily/atmos_daily_together_2_5.nc', chunks={'time': 30});     data_vdfmr = data_vdfmr['var249'].load().loc['1979-01':'2016-12'].sel(plev=85000.)
#ds['vdfmr'] = pentad_mean_climatology(data_vdfmr, np.arange(1979,2017))
#data_admr = xr.open_dataset('/disca/share/rg419/JRA_55/admr_daily/atmos_daily_together_2_5.nc', chunks={'time': 30});     data_admr = data_admr['var236'].load().loc['1979-01':'2016-12'].sel(plev=85000.)
#ds['admr'] = pentad_mean_climatology(data_admr, np.arange(1979,2017))

#print(ds)

#data_t = xr.open_dataset('/disca/share/rg419/jra_temp_daily_850.nc', chunks={'time': 30});     data_t = data_t['var11'].load().loc['1979-01':'2016-12'].sel(lev=85000.)
#data_q = xr.open_dataset('/disca/share/rg419/jra_sphum_daily_850.nc', chunks={'time': 30});     data_q = data_q['var51'].load().loc['1979-01':'2016-12'].sel(lev=85000.)
#data_h = xr.open_dataset('/disca/share/rg419/jra_height_daily_850.nc', chunks={'time': 30});     data_h = data_h['var7'].load().loc['1979-01':'2016-12'].sel(lev=85000.)

data_u = xr.open_dataset('/disca/share/rg419/jra_ucomp_daily_850.nc', chunks={'time': 30});     data_u = data_u['var33'].load().loc['1979-01':'2016-12']
data_w = xr.open_dataset('/disca/share/rg419/jra_omega_daily_850.nc', chunks={'time': 30});     data_w = data_w['var39'].load().loc['1979-01':'2016-12'].sel(lev=85000.)
data_v_temp = xr.open_dataset('/disca/share/rg419/jra_vcomp_daily_850.nc', chunks={'time': 30});    data_v_temp = data_v_temp['var34'].load().loc['1979-01':'2016-12']
data_v = xr.DataArray(data_v_temp.sel(lev=85000.).values, coords={'time': data_u.time, 'lat': data_u.lat, 'lon': data_u.lon}, dims=('time','lat','lon'))


#ds['ucomp_temp'] = pentad_mean_climatology( xr.DataArray(data_t.values * data_u.values, coords={'time': data_u.time, 'lat': data_u.lat, 'lon': data_u.lon}, dims=('time','lat','lon'))  , np.arange(1979,2017))
#ds['vcomp_temp'] = pentad_mean_climatology( xr.DataArray(data_t.values * data_v.values, coords={'time': data_u.time, 'lat': data_u.lat, 'lon': data_u.lon}, dims=('time','lat','lon'))  , np.arange(1979,2017))
#ds['omega_temp'] = pentad_mean_climatology( xr.DataArray(data_t.values * data_w.values, coords={'time': data_u.time, 'lat': data_u.lat, 'lon': data_u.lon}, dims=('time','lat','lon'))  , np.arange(1979,2017))

#ds['ucomp_height'] = pentad_mean_climatology( xr.DataArray(data_h.values * data_u.values, coords={'time': data_u.time, 'lat': data_u.lat, 'lon': data_u.lon}, dims=('time','lat','lon'))  , np.arange(1979,2017))
#ds['vcomp_height'] = pentad_mean_climatology( xr.DataArray(data_h.values * data_v.values, coords={'time': data_u.time, 'lat': data_u.lat, 'lon': data_u.lon}, dims=('time','lat','lon'))  , np.arange(1979,2017))
#ds['omega_height'] = pentad_mean_climatology( xr.DataArray(data_h.values * data_w.values, coords={'time': data_u.time, 'lat': data_u.lat, 'lon': data_u.lon}, dims=('time','lat','lon'))  , np.arange(1979,2017))

#ds['sphum_u'] = pentad_mean_climatology( xr.DataArray(data_q.values * data_u.values, coords={'time': data_u.time, 'lat': data_u.lat, 'lon': data_u.lon}, dims=('time','lat','lon'))  , np.arange(1979,2017))
#ds['sphum_v'] = pentad_mean_climatology( xr.DataArray(data_q.values * data_v.values, coords={'time': data_u.time, 'lat': data_u.lat, 'lon': data_u.lon}, dims=('time','lat','lon'))  , np.arange(1979,2017))
#ds['sphum_w'] = pentad_mean_climatology( xr.DataArray(data_q.values * data_w.values, coords={'time': data_u.time, 'lat': data_u.lat, 'lon': data_u.lon}, dims=('time','lat','lon'))  , np.arange(1979,2017))

ds['ucomp'] = pentad_mean_climatology( data_u, np.arange(1979,2017))
ds['vcomp'] = pentad_mean_climatology( data_v, np.arange(1979,2017))
ds['omega'] = pentad_mean_climatology( data_w, np.arange(1979,2017))

#ds['temp'] = pentad_mean_climatology(data_t, np.arange(1979,2017))
#ds['sphum'] = pentad_mean_climatology(data_q, np.arange(1979,2017))
#ds['height'] = pentad_mean_climatology(data_h, np.arange(1979,2017))

print(ds)

fileout = '/disca/share/rg419/JRA_55/heating_terms3.nc'  
ds.to_netcdf(path=fileout)

