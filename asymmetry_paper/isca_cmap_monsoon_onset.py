''' 
5/03/2019 Make onset date figure: CMAP in NH (and SH upside down?) Isca: half_land, half_land_nh_tibet
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


def get_onset(rain_rel): 
    rain_rel_masked = np.ma.masked_less(rain_rel.values, 5)  # Mask relative rainfall where it is less than 5mm/day
    onset_index = np.ma.notmasked_edges(rain_rel_masked, axis=0) # Look along pentad axis to find edges of mask
    onset = np.zeros((len(rain_rel.lat),len(rain_rel.lon))) # Create array of nans to load onsets into
    onset[:] = np.nan
    
    for i in range(0,len(onset_index[0][1])): # Extract onsets from mask edges
        onset[ onset_index[0][1][i], onset_index[0][2][i] ] = onset_index[0][0][i]+1
    onset_pentad = xr.DataArray(onset, [('lat', rain_rel.lat), ('lon', rain_rel.lon)]) # Make onset pentad into a dataarray
    return onset_pentad
    

def relative_rain_isca(run, ax, land_mask_name=None):
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    jan_precip = data.precipitation.sel(xofyear=range(1,7)).mean('xofyear')
    relative_rain = (data.precipitation - jan_precip)*86400.
    onset_pentad = get_onset(relative_rain)
        
    f1 = onset_pentad.plot.pcolormesh(ax=ax, x='lon',y='lat',  levels=np.arange(25.,56.), cmap='plasma_r', add_colorbar=False, add_labels=False)
        
    if land_mask_name==None:
        land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/' + run + '.nc'
    else:
        land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/' + land_mask_name + '.nc'
    
    land = xr.open_dataset(land_mask)
    land.land_mask.plot.contour(ax=ax, x='lon', y='lat', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
    land.zsurf.plot.contour(ax=ax, x='lon', y='lat', levels=np.arange(2000.,3001.,1000.), add_labels=False, colors='0.6')
    ax.grid(True,linestyle=':')
    ax.set_ylim([0,45])
    ax.set_xticks(np.arange(0.,361.,60.))
    ax.set_yticks(np.arange(0.,46.,15.))
    ax.set_ylabel('Latitude')
    
    data.close()
    return f1


def relative_rain_cmap(ax, nh=True):
    '''10/10/2018 Load up CMAP pentad and monthly means, and use to evaluate Northern and Southern hemisphere monsoon onset using the Wang and LinHo criteria'''
    
    # Load up data
    data = xr.open_dataset('/disca/share/rg419/CMAP_precip.pentad.mean.nc')
    data_mon = xr.open_dataset('/disca/share/rg419/CMAP_precip.mon.mean.nc')
    
    # Add pentad as a coordinate and make a climatology
    data.coords['pentad'] = (('time'), np.tile(np.arange(1,74),38))
    precip_clim_mon = data_mon.precip.groupby('time.month').mean('time')
    precip_clim_pentad = data.precip.groupby('pentad').mean('time')
    
    # Calculate rainfal relative to January and July. For July re-order pentads so onset can be identified for Southern hemisphere monsoons too.
    if nh:
        rain_rel = xr.DataArray((precip_clim_pentad - precip_clim_mon.sel(month=1)).values, [('pentad', np.arange(1,74)), ('lat', data.lat), ('lon', data.lon)])
        onset = get_onset(rain_rel)
        f1 = onset.plot.pcolormesh(ax=ax, x='lon',y='lat',  levels=np.arange(20.,51.), cmap='plasma_r', add_colorbar=False, add_labels=False)
        ax.set_ylim([0,45])
        ax.set_yticks(np.arange(0.,46.,15.))   
    else:
        rain_rel = xr.DataArray((precip_clim_pentad - precip_clim_mon.sel(month=7)).values, [('pentad', (np.arange(1,74)-37)%73-36), ('lat', data.lat), ('lon', data.lon)])
        rain_rel = rain_rel.sortby('pentad')        
        onset= get_onset(rain_rel) - 37.   # Onset is output as an index along the time axis. For July, where this has been rotated, therefore subtract 37 to correct
        f1 = onset.plot.pcolormesh(ax=ax, x='lon',y='lat', levels=np.arange(-15.,16.), cmap='plasma_r', add_colorbar=False, add_labels=False, yincrease=False)
        ax.set_ylim([-45,0])
        ax.set_yticks(np.arange(-45.,1.,15.))
    
    land_mask = '/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc'
    land = xr.open_dataset(land_mask)
    land.lsm[0,:,:].plot.contour(ax=ax, x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
    (land.z[0,:,:]/9.81).plot.contour(ax=ax, x='longitude', y='latitude', levels=np.arange(2000.,3001.,1000.), add_labels=False, colors='0.6')
    
    ax.grid(True,linestyle=':')
    ax.set_ylabel('Latitude')    
    ax.set_xticks(np.arange(0.,361.,60.))
    
    data.close()
    data_mon.close()
    return f1



# Set figure parameters
rcParams['figure.figsize'] = 7, 5
rcParams['font.size'] = 10

# Start figure with 4 subplots
fig, ((ax1, ax2, ax3, ax4)) = plt.subplots(4, 1, sharex='col')
axes = [ax1, ax2, ax3, ax4]

f1 = relative_rain_isca('half_land', ax1, land_mask_name='half_shallow')
f2 = relative_rain_isca('half_nh_land_tibet_ol8', ax2, land_mask_name='half_nh_land_tibet')
f3 = relative_rain_cmap(ax3, nh=True)
f4 = relative_rain_cmap(ax4, nh=False)

ax4.set_xlabel('Longitude')

plt.subplots_adjust(left=0.12, right=0.92, top=0.97, bottom=0.1, hspace=0.2)
cb1=fig.colorbar(f1, ax=ax1, ticks=np.arange(25.,56.,10.), use_gridspec=True, orientation = 'vertical',fraction=0.05, pad=0.03, aspect=15, shrink=1)
#cb1.set_label('Onset Pentad')
cb1=fig.colorbar(f2, ax=ax2, ticks=np.arange(25.,56.,10.), use_gridspec=True, orientation = 'vertical',fraction=0.05, pad=0.03, aspect=15, shrink=1)
#cb1.set_label('Onset Pentad')
cb1=fig.colorbar(f3, ax=ax3, ticks=np.arange(20.,51.,10.), use_gridspec=True, orientation = 'vertical',fraction=0.05, pad=0.03, aspect=15, shrink=1)
#cb1.set_label('Onset Pentad')
cb1=fig.colorbar(f4, ax=ax4, ticks=np.arange(-15.,16.,10.), use_gridspec=True, orientation = 'vertical',fraction=0.05, pad=0.03, aspect=15, shrink=1)
cb1.ax.set_yticklabels((np.arange(-15,16,10)+73)%73)  # Label colorbar with pentad numbers corresponding to the onset pentad
#cb1.set_label('Onset Pentad')

ax4.invert_yaxis()

# Save as a pdf
plt.savefig(plot_dir + 'monsoon_onset_isca_cmap.pdf', format='pdf')
plt.close()

