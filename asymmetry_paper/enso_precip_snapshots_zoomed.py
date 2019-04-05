''' 
05/04/2019 Plot up snapshots of precip and lower level wind for climatology, El Nino, Neutral, and La Nina years
'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
import pandas as pd

# Set plotting directory
plot_dir = '/scratch/rg419/plots/asymmetry_paper/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


def pentad_mean(data, years):  # Function to convert to pentad means
    # Define empty array to put pentads into
    pentad_years = np.array([])
    i=0
    for year in years:
        #Select daily datapoints in a given year
        data_year = data.sel(time=str(year))
        # Count these. If a leap year add an extra day to pentad 12 (to account for 29th Feb)
        if len(data_year.time)==366:
            pentad = np.repeat(np.arange(1., 74.), 5)
            pentad = np.insert(pentad, 59, 12)   
            print(pentad) 
        else:
            pentad = np.repeat(np.arange(1., 74.), 5)    
        # To resample as pentad means, add 100*yearno to distinguish years
        pentad_years = np.concatenate((pentad_years, pentad + i*100.))
        i=i+1
    
    # Add pentad coordinate to data    
    data = data.assign_coords(pentad = ('time', pentad_years))
    
    # Groupby pentad and take time mean
    data_pentads = data.groupby('pentad').mean(('time'))
    # Now remove the 100*yearno factor
    data_pentads['pentad'] = data_pentads['pentad'] - np.repeat(np.arange(0.,3800.,100.), 73)

    return data_pentads
    
    
# Decide whether nino index is El Nino or La Nina
def get_phase(nino_34):
    enso_phase = []
    i=0
    for i in range(len(nino_34.time)):
        if nino_34[i] >= 0.5:
            enso_phase.append('El Nino')
        elif nino_34[i] <= -0.5:
            enso_phase.append('La Nina')
        else:
            enso_phase.append('Neutral')
        i=i+1    
    nino_34_out = xr.DataArray(nino_34, coords={'enso_phase': ('time', enso_phase), 'time': ('time', nino_34.time)}, dims=['time'])
    return nino_34_out


def plot_gill_dev_cmap(ax, pentad):
    
    nino_34 = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/nino_34.nc')
    nino_34_march = nino_34.nino_34.sel(time=pd.to_datetime(['%04d-03-15' % y for y in range(1979,2017)]))
    nino_34_march = get_phase(nino_34_march)
    
    data_precip = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/CMAP_precip.pentad.mean.nc', chunks={'time': 30})
    data_precip = data_precip.load()
    data_slp = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/jra_slp_daily.nc', chunks={'time': 30})
    data_slp = data_slp.load()
    data_slp = data_slp.sel(time=slice('1979','2016'))
    data_slp = pentad_mean(data_slp.var2/100., np.arange(1979,2017))
    a=boobs
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
#pentads=[20,26,28,32,40,47]
#pentads=[20,26,32,38,44,50]
pentads=[22,28,34,40,46]
for ax in axes[::3]:
    f1 = plot_gill_dev_cmap(ax, pentads[i])
    i=i+1

for ax in [ax13, ax14, ax15]:
    ax.set_xlabel('Longitude')

plt.subplots_adjust(left=0.07, right=0.97, top=0.97, bottom=0.05, hspace=0.3, wspace=0.15)
cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.02, pad=0.07, aspect=30, shrink=0.5)
cb1.set_label('Precipitation, mm/day')

# Save as a pdf
plt.savefig(plot_dir + 'enso_precip_snapshots_zoomed.pdf', format='pdf')
plt.close()

