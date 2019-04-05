'''6/12/2018 Plot pentad by pentad development of JRA Rossby no and precip
'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from windspharm.xarray import VectorWind
import matplotlib.patches as patches
from data_handling_updates import gradients as gr, model_constants as mc

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

data = xr.open_dataset('/disca/share/sit204/jra_55/1958_2016/omega_daily/atmos_daily_together.nc', chunks={'time': 30})
data = pentad_mean_climatology(data, np.arange(1958,2017))
data.to_netcdf('/disca/share/rg419/jra_omega_pentad_clim_alllevs.nc')
#data.to_netcdf('/disca/share/rg419/jra_omega_daily_200.nc')


#data_u = xr.open_dataset('/disca/share/sit204/jra_55/1958_2016/ucomp_daily/atmos_daily_together.nc', chunks={'time': 30});
#data_u = pentad_mean_climatology(data_u, np.arange(1958,2017))
#data_u.to_netcdf('/disca/share/rg419/jra_ucomp_pentad_clim.nc')

#data_w = xr.open_dataset('/disca/share/rg419/jra_omega_daily_200.nc', chunks={'time': 30});
#data_w = pentad_mean_climatology(data_w, np.arange(1958,2017))
#data_w.to_netcdf('/disca/share/rg419/jra_omega_pentad_clim.nc')

#data_v = xr.open_dataset('/disca/share/rg419/jra_vcomp_daily_200.nc', chunks={'time': 30})
#data_v = pentad_mean_climatology(data_v, np.arange(1958,2017))
#data_v.to_netcdf('/disca/share/rg419/jra_vcomp_pentad_clim.nc')