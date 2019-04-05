''' 
12/02/2019 Make comparable hms: half_land, half_10_land, half_10_land_tibet, and CMAP
'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from data_handling_updates import model_constants as mc, gradients as gr

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
    data['mse'] = (mc.cp_air * data.temp + mc.L * data.sphum + mc.grav * data.height)/1000.
    data['mse_u'] = (mc.cp_air * data.ucomp_temp + mc.L * data.sphum_u + mc.grav * data.ucomp_height)/1000.
    data['mse_v'] = (mc.cp_air * data.vcomp_temp + mc.L * data.sphum_v + mc.grav * data.vcomp_height)/1000.
    
    data['mse_v_div'] = -1.*(gr.ddx(data.mse_u) + gr.ddy(data.mse_v))
    
    # Take zonal anomaly
    data_zanom = data - data.mean('lon')
    
    title = 'Pentad ' + str(int(pentad))
    
    f1 = data.mse_v_div.sel(xofyear=pentad, pfull=850.).plot.contourf(x='lon', y='lat', ax=ax, add_labels=False, add_colorbar=False, extend='both', zorder=1, levels = np.arange(-0.004, 0.0041, 0.0005))
    #f1 = data.mse.sel(xofyear=pentad, pfull=850.).plot.contourf(x='lon', y='lat', ax=ax, add_labels=False, add_colorbar=False, extend='both',  cmap='Blues', zorder=1, levels = np.arange(310.,351.,2.5))
    ax.contour(data.lon, data.lat, data.mse.sel(pfull=850.,xofyear=pentad), levels = np.arange(310.,351.,2.5), colors='0.4', zorder=2)
    #ax.contour(data.lon, data.lat, data.mse_v_div.sel(pfull=850.,xofyear=pentad), levels = np.arange(0., 0.0041, 0.0005), colors='0.4', zorder=2)
    #ax.contour(data.lon, data.lat, data.mse_v_div.sel(pfull=850.,xofyear=pentad), levels = np.arange(-0.004, 0., 0.0005), colors='0.4', linestyles='--', zorder=2)
    b = ax.quiver(data.lon[::6], data.lat[::3], data.mse_u.sel(pfull=850.,xofyear=pentad)[::3,::6], data.mse_v.sel(pfull=850.,xofyear=pentad)[::3,::6], scale=50000, angles='xy', width=0.01, headwidth=3., headlength=5., zorder=3)
    #ax.quiverkey(b, 90.,35., 5, str(5) + ' m/s', coordinates='data', fontproperties={'weight': 'bold', 'size': 10}, color='k', labelcolor='k', labelsep=0.03, zorder=10)
    ax.grid(True,linestyle=':')
    ax.set_ylim(-15.,45.)
    ax.set_yticks(np.arange(-15.,45.,15.))
    ax.set_xlim(90.,225.)
    ax.set_xticks(np.arange(90.,226.,45.))
    ax.set_title(title, fontsize=11)
    land = xr.open_dataset(land_mask)
    land.land_mask.plot.contour(x='lon', y='lat', ax=ax, levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')    
    land.zsurf.plot.contour(ax=ax, x='lon', y='lat', levels=np.arange(2000.,3001.,1000.), add_labels=False, colors='k')
    
    return f1, b



def plot_gill_dev_cmap(ax, pentad):
   
    data_t = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/jra_temp_daily_850.nc', chunks={'time': 30});     data_t = data_t['var11'].load().loc['1979-01':'2016-12']
    data_q = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/jra_sphum_daily_850.nc', chunks={'time': 30});     data_q = data_q['var51'].load().loc['1979-01':'2016-12']
    data_h = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/jra_height_daily_850.nc', chunks={'time': 30});     data_h = data_h['var7'].load().loc['1979-01':'2016-12']

    data_u = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/jra_ucomp_daily_850.nc', chunks={'time': 30});     data_u = data_u['var33'].load().loc['1979-01':'2016-12']
    data_v_temp = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/jra_vcomp_daily_850.nc', chunks={'time': 30});    data_v_temp = data_v_temp['var34'].load().loc['1979-01':'2016-12']
    data_v = xr.DataArray(data_v_temp.sel(lev=85000.).values, coords={'time': data_u.time, 'lat': data_u.lat, 'lon': data_u.lon}, dims=('time','lat','lon'))    
    
    data_mse = (mc.cp_air * data_t + mc.L * data_q + 9.81 * data_h)/1000.
    data_mse_u = xr.DataArray(data_mse.sel(lev=85000.).values * data_u.values, coords={'time': data_u.time, 'lat': data_u.lat, 'lon': data_u.lon}, dims=('time','lat','lon'))    
    data_mse_v = xr.DataArray(data_mse.sel(lev=85000.).values * data_v.values, coords={'time': data_u.time, 'lat': data_u.lat, 'lon': data_u.lon}, dims=('time','lat','lon'))    

    data_u = pentad_mean_climatology(data_u, np.arange(1979,2017))
    data_v = pentad_mean_climatology(data_v, np.arange(1979,2017))
    data_t = pentad_mean_climatology(data_t, np.arange(1979,2017))
    data_q = pentad_mean_climatology(data_q, np.arange(1979,2017))
    data_h = pentad_mean_climatology(data_h, np.arange(1979,2017))
    data_mse = pentad_mean_climatology(data_mse, np.arange(1979,2017))
    data_mse_u = pentad_mean_climatology(data_mse_u, np.arange(1979,2017))
    data_mse_v = pentad_mean_climatology(data_mse_v, np.arange(1979,2017))
    
    data_mse_v_div = -1.*(gr.ddx(data_mse_u) + gr.ddy(data_mse_v))
    data_u = data_u - data_u.mean('lon')
    data_v = data_v - data_v.mean('lon')
    #data_t = data_t - data_t.mean('lon')
    #data_q = data_q - data_q.mean('lon')
    #data_h = data_h - data_h.mean('lon')
    title = 'Pentad ' + str(int(pentad))      
    f1 = data_mse_v_div.sel(pentad=pentad).plot.contourf(x='lon', y='lat', ax=ax, add_labels=False, add_colorbar=False, extend='both', zorder=1, levels = np.arange(-0.004, 0.0041, 0.0005))
    #f1 = data_mse.sel(pentad=pentad,lev=85000.).plot.contourf(x='lon', y='lat', ax=ax, add_labels=False, add_colorbar=False, extend='both',  cmap='Blues', zorder=1, levels = np.arange(310.,351.,2.5))
    ax.contour(data_mse.lon, data_mse.lat, data_mse.sel(pentad=pentad,lev=85000.), levels = np.arange(310.,351.,2.5), colors='0.4', zorder=2)
    #ax.contour(data_mse.lon, data_mse.lat, data_mse_v_div.sel(pentad=pentad), levels = np.arange(0., 0.0041, 0.0005), colors='0.4', zorder=2)
    #ax.contour(data_mse.lon, data_mse.lat, data_mse_v_div.sel(pentad=pentad), levels = np.arange(-0.004, 0., 0.0005), colors='0.4', linestyles='--', zorder=2)
    b = ax.quiver(data_u.lon[::6], data_u.lat[::3], data_mse_u.sel(pentad=pentad)[::3,::6], data_mse_v.sel(pentad=pentad)[::3,::6], scale=50000, angles='xy', width=0.01, headwidth=3., headlength=5., zorder=3)
    #ax.quiverkey(b, 90.,35, 5, str(5) + ' m/s', fontproperties={'weight': 'bold', 'size': 10}, color='k', labelcolor='k', labelsep=0.03, zorder=10)
    
    ax.grid(True,linestyle=':')
    ax.set_ylim(-15.,45.)
    ax.set_yticks(np.arange(-15.,45.,15.))
    ax.set_xlim(45.,180.)
    ax.set_xticks(np.arange(45.,181.,45.))
    ax.set_title(title, fontsize=11)
    land_mask = '/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc'
    land = xr.open_dataset(land_mask)
    land.lsm[0,:,:].plot.contour(ax=ax, x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
    (land.z[0,:,:]/9.81).plot.contour(ax=ax, x='longitude', y='latitude', levels=np.arange(2000.,3001.,1000.), add_labels=False, colors='k')
    
    return f1, b


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
    f1, b = plot_gill_dev_isca('half_land', ax, pentads[i], land_mask='/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')
    ax.set_ylabel('Latitude')
    i=i+1

i=0
for ax in axes[1::3]:
    f1, b = plot_gill_dev_isca('half_nh_land_tibet_ol8', ax, pentads[i], land_mask='/scratch/rg419/Experiments/asym_aquaplanets/input/half_nh_land_tibet.nc')
    i=i+1

i=0
#pentads=[20,26,28,32,40,47]
#pentads=[20,26,32,38,44,50]
pentads=[22,28,34,40,46]
for ax in axes[2::3]:
    f1, b = plot_gill_dev_cmap(ax, pentads[i])
    i=i+1

for ax in [ax13, ax14, ax15]:
    ax.set_xlabel('Longitude')

ax1.quiverkey(b, 135.,0., 5., str(5) + ' m/s', fontproperties={'weight': 'bold', 'size': 10}, color='k', labelcolor='k', labelsep=0.03, zorder=10)

plt.subplots_adjust(left=0.07, right=0.97, top=0.97, bottom=0.05, hspace=0.3, wspace=0.15)
cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.02, pad=0.07, aspect=30, shrink=0.5)
cb1.set_label('Moist Static Energy, kJ/kg')

# Save as a pdf
plt.savefig(plot_dir + 'mseflux_snapshots_zoomed.pdf', format='pdf')
plt.close()

