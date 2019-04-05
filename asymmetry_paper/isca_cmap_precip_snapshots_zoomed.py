''' 
12/02/2019 Make comparable hms: half_land, half_10_land, half_10_land_tibet, and CMAP
'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh

# Set plotting directory
plot_dir = '/scratch/rg419/plots/asymmetry_paper/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


def pentad_mean_climatology(data, years):  # Function to get pentad of year
    pentad_years = np.array([])
    for year in years:
        data_year = data.sel(time=str(year))
        if len(data_year.time)==366:
            pentad = np.repeat(np.arange(1., 74.), 5)
            pentad = np.insert(pentad, 59, 12)    
        else:
            pentad = np.repeat(np.arange(1., 74.), 5)    
        pentad_years = np.concatenate((pentad_years, pentad))
        
    data = data.assign_coords(pentad = ('time', pentad_years))
    
    data_pentads = data.groupby('pentad').mean(('time'))
    
    return data_pentads


def plot_gill_dev_isca(run, ax, pentad, land_mask=None):
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    data['precipitation'] = (data.precipitation*86400.)
    # Take zonal anomaly
    data_zanom = data - data.mean('lon')
    
    title = 'Pentad ' + str(int(pentad))
    f1 = data.precipitation.sel(xofyear=pentad).plot.contourf(x='lon', y='lat', ax=ax, levels = np.arange(3.,19.,3.), add_labels=False, add_colorbar=False, extend='both',  cmap='Blues', zorder=1)
    ax.contour(data_zanom.lon, data_zanom.lat, data_zanom.slp.sel(xofyear=pentad), levels = np.arange(0.,16.,3.), colors='0.4', alpha=0.5, zorder=3)
    ax.contour(data_zanom.lon, data_zanom.lat, data_zanom.slp.sel(xofyear=pentad), levels = np.arange(-15.,0.,3.), colors='0.4', alpha=0.5, linestyle='--', zorder=3)
    b = ax.quiver(data.lon[::6], data.lat[::3], data_zanom.ucomp.sel(pfull=850.,xofyear=pentad)[::3,::6], data_zanom.vcomp.sel(pfull=850.,xofyear=pentad)[::3,::6], scale=100, angles='xy', width=0.01, headwidth=3., headlength=5., zorder=3)
    ax.grid(True,linestyle=':')
    ax.set_ylim(-15.,45.)
    ax.set_yticks(np.arange(-15.,45.,15.))
    ax.set_xlim(90.,225.)
    ax.set_xticks(np.arange(90.,226.,45.))
    ax.set_title(title, fontsize=11)
    land = xr.open_dataset(land_mask)
    land.land_mask.plot.contour(x='lon', y='lat', ax=ax, levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', zorder=2)    
    land.zsurf.plot.contour(ax=ax, x='lon', y='lat', levels=np.arange(2000.,3001.,1000.), add_labels=False, colors='k', zorder=2)
    
    return f1



def plot_gill_dev_cmap(ax, pentad):
    
    data_slp = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/jra_slp_daily.nc', chunks={'time': 30})
    data_precip = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/CMAP_precip.pentad.mean.nc', chunks={'time': 30})
    data_precip = data_precip.load()
    data_slp = data_slp.load()
    data_slp = data_slp.sel(time=slice('1979','2016'))
    data_slp = pentad_mean_climatology(data_slp.var2/100., np.arange(1979,2017))
    data_precip.coords['pentad'] = (('time'), np.tile(np.arange(1,74),38))
    data_precip = data_precip.groupby('pentad').mean('time')
    data_slp = data_slp - data_slp.mean('lon')
    
    data_u = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/jra_ucomp_daily_850.nc', chunks={'time': 30});     data_u = data_u['var33'].load().loc['1979-01':'2016-12']
    data_v_temp = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/jra_vcomp_daily_850.nc', chunks={'time': 30});    data_v_temp = data_v_temp['var34'].load().loc['1979-01':'2016-12']
    data_v = xr.DataArray(data_v_temp.sel(lev=85000.).values, coords={'time': data_u.time, 'lat': data_u.lat, 'lon': data_u.lon}, dims=('time','lat','lon'))    
    data_u = pentad_mean_climatology(data_u, np.arange(1979,2017))
    data_v = pentad_mean_climatology(data_v, np.arange(1979,2017))
    data_u = data_u - data_u.mean('lon')
    data_v = data_v - data_v.mean('lon')
    
    title = 'Pentad ' + str(int(pentad))       
    f1 = data_precip.precip.sel(pentad=pentad).plot.contourf(x='lon', y='lat', ax=ax, levels = np.arange(3.,19.,3.), add_labels=False, add_colorbar=False, extend='both',  cmap='Blues', zorder=1)
    ax.contour(data_slp.lon, data_slp.lat, data_slp.sel(pentad=pentad), levels = np.arange(0.,16.,3.), colors='0.4', alpha=0.5, zorder=3)
    ax.contour(data_slp.lon, data_slp.lat, data_slp.sel(pentad=pentad), levels = np.arange(-15.,0.,3.), colors='0.4', alpha=0.5, linestyle='--', zorder=3)
    b = ax.quiver(data_u.lon[::6], data_u.lat[::3], data_u.sel(pentad=pentad)[::3,::6], data_v.sel(pentad=pentad)[::3,::6], scale=100, angles='xy', width=0.01, headwidth=3., headlength=5., zorder=3)
    ax.grid(True,linestyle=':')
    ax.set_ylim(-15.,45.)
    ax.set_yticks(np.arange(-15.,45.,15.))
    ax.set_xlim(45.,180.)
    ax.set_xticks(np.arange(45.,181.,45.))
    ax.set_title(title, fontsize=11)
    land_mask = '/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc'
    land = xr.open_dataset(land_mask)
    land.lsm[0,:,:].plot.contour(ax=ax, x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', zorder=2)
    (land.z[0,:,:]/9.81).plot.contour(ax=ax, x='longitude', y='latitude', levels=np.arange(2000.,3001.,1000.), add_labels=False, colors='k', zorder=2)
    
    return f1


# Set figure parameters
rcParams['figure.figsize'] = 10, 11
rcParams['font.size'] = 11

# Start figure with 18 subplots
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9), (ax10, ax11, ax12), (ax13, ax14, ax15)) = plt.subplots(5, 3, sharey='row', sharex='col')
axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12, ax13, ax14, ax15]

i=0
#pentads=[30,36,38,42,50,57]
#pentads=[30,36,42,48,54,60]
pentads=[32,38,44,50,56]
for ax in axes[::3]:
    f1 = plot_gill_dev_isca('half_land', ax, pentads[i], land_mask='/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')
    ax.set_ylabel('Latitude')
    i=i+1

i=0
for ax in axes[1::3]:
    f1 = plot_gill_dev_isca('half_nh_land_tibet_ol8', ax, pentads[i], land_mask='/scratch/rg419/Experiments/asym_aquaplanets/input/half_nh_land_tibet.nc')
    i=i+1

i=0
#pentads=[20,26,28,32,40,47]
#pentads=[20,26,32,38,44,50]
pentads=[22,28,34,40,46]
for ax in axes[2::3]:
    f1 = plot_gill_dev_cmap(ax, pentads[i])
    i=i+1

for ax in [ax13, ax14, ax15]:
    ax.set_xlabel('Longitude')

plt.subplots_adjust(left=0.07, right=0.97, top=0.97, bottom=0.05, hspace=0.3, wspace=0.15)
cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.02, pad=0.07, aspect=30, shrink=0.5)
cb1.set_label('Precipitation, mm/day')

# Save as a pdf
plt.savefig(plot_dir + 'precip_snapshots_zoomed.pdf', format='pdf')
plt.close()

