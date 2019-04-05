'''5/12/2018 As gill_development.py but for data from JRA-55 and CMAP precip
'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from windspharm.xarray import VectorWind
from data_handling_updates import model_constants as mc

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
    
def plot_gill_dev(lev=85000, qscale=100., windtype='full', ref_arrow=5, mse=False, video=False):
    
    if mse:
        data_t = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/jra_temp_daily_850.nc', chunks={'time': 30});     data_t = data_t['var11'].load().loc['1979-01':'2016-12']
        data_q = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/jra_sphum_daily_850.nc', chunks={'time': 30});     data_q = data_q['var51'].load().loc['1979-01':'2016-12']
        data_h = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/jra_height_daily_850.nc', chunks={'time': 30});     data_h = data_h['var7'].load().loc['1979-01':'2016-12']
        data_t = pentad_mean_climatology(data_t, np.arange(1979,2017))
        data_q = pentad_mean_climatology(data_q, np.arange(1979,2017))
        data_h = pentad_mean_climatology(data_h, np.arange(1979,2017))
        data_mse = (mc.cp_air * data_t + mc.L * data_q + 9.81 * data_h)/1000.
    else:
        data_precip = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/CMAP_precip.pentad.mean.nc', chunks={'time': 30})
        data_precip = data_precip.load()
        data_precip.coords['pentad'] = (('time'), np.tile(np.arange(1,74),38))
        data_precip = data_precip.groupby('pentad').mean('time')
    
    if windtype is not 'none':    
        data_u = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/jra_ucomp_daily_850.nc', chunks={'time': 30});     data_u = data_u['var33'].load().loc['1979-01':'2016-12']
        data_v_temp = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/jra_vcomp_daily_850.nc', chunks={'time': 30});    data_v_temp = data_v_temp['var34'].load().loc['1979-01':'2016-12']
        # v has different time coord to u, presumably due to how Stephen has downloaded/averaged. I think the two are equivalent, so just substitute the time dimension into v
        data_v = xr.DataArray(data_v_temp.sel(lev=85000.).values, coords={'time': data_u.time, 'lat': data_u.lat, 'lon': data_u.lon}, dims=('time','lat','lon'))
        data_u = pentad_mean_climatology(data_u, np.arange(1979,2017))
        data_v = pentad_mean_climatology(data_v, np.arange(1979,2017))
        data_u = data_u - data_u.mean('lon')
        data_v = data_v - data_v.mean('lon')
        if windtype is not 'full':
            # Get rotational and divergent components of the flow
            w = VectorWind(data_u, data_v)
            streamfun, vel_pot = w.sfvp()
            uchi, vchi, upsi, vpsi = w.helmholtz()
            #print(uchi.lat)
            #print(data_u.lat)
            uchi_zanom = (uchi - uchi.mean('lon')).sortby('lat')
            vchi_zanom = (vchi - vchi.mean('lon')).sortby('lat')
            upsi_zanom = (upsi - upsi.mean('lon')).sortby('lat')
            vpsi_zanom = (vpsi - vpsi.mean('lon')).sortby('lat')
        
    data_slp = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/jra_slp_daily.nc', chunks={'time': 30});    data_slp = data_slp['var2'].load().loc['1979-01':'2016-12']
    data_slp = data_slp.load()
    data_slp = pentad_mean_climatology(data_slp/100., np.arange(1979,2017))
    data_slp = data_slp - data_slp.mean('lon')
    print('files opened')
    
    
    rcParams['figure.figsize'] = 10, 5
    rcParams['font.size'] = 14
    
    for i in range(73):        
        fig, ax1 = plt.subplots()
        title = 'Pentad ' + str(int(data_u.pentad[i]))
        
        if mse:
            f1 = data_mse.sel(pentad=i+1, lev=85000.).plot.contourf(x='lon', y='lat', ax=ax1, add_labels=False, add_colorbar=False, extend='both',  cmap='Blues', zorder=1, levels = np.arange(290.,341.,5.))
        else:
            f1 = data_precip.precip[i,:,:].plot.contourf(x='lon', y='lat', ax=ax1, levels = np.arange(2.,15.,2.), add_labels=False, add_colorbar=False, extend='both',  cmap='Blues', zorder=1)

        ax1.contour(data_slp.lon, data_slp.lat, data_slp[i,:,:], levels = np.arange(0.,16.,3.), colors='0.4', alpha=0.5, zorder=2)
        ax1.contour(data_slp.lon, data_slp.lat, data_slp[i,:,:], levels = np.arange(-15.,0.,3.), colors='0.4', alpha=0.5, linestyle='--', zorder=2)
        if windtype=='div':
            b = ax1.quiver(uchi_zanom.lon[::6], uchi_zanom.lat[::3], uchi_zanom[i,::3,::6], vchi_zanom[i,::3,::6], scale=qscale, angles='xy', width=0.005, headwidth=3., headlength=5., zorder=3)
            ax1.quiverkey(b, 5.,65., ref_arrow, str(ref_arrow) + ' m/s', fontproperties={'weight': 'bold', 'size': 10}, color='k', labelcolor='k', labelsep=0.03, zorder=10)
        elif windtype=='rot':
            b = ax1.quiver(upsi_zanom.lon[::6], upsi_zanom.lat[::3], upsi_zanom[i,::3,::6], vpsi_zanom[i,::3,::6], scale=qscale, angles='xy', width=0.005, headwidth=3., headlength=5., zorder=3)
            ax1.quiverkey(b, 5.,65., ref_arrow, str(ref_arrow) + ' m/s', fontproperties={'weight': 'bold', 'size': 10}, color='k', labelcolor='k', labelsep=0.03, zorder=10)
        elif windtype=='full':
            b = ax1.quiver(data_u.lon[::6], data_u.lat[::3], data_u[i,::3,::6], data_v[i,::3,::6], scale=qscale, angles='xy', width=0.005, headwidth=3., headlength=5., zorder=3)
            ax1.quiverkey(b, 5.,65., ref_arrow, str(ref_arrow) + ' m/s', fontproperties={'weight': 'bold', 'size': 10}, color='k', labelcolor='k', labelsep=0.03, zorder=10)
        else:
            windtype='none'
        ax1.grid(True,linestyle=':')
        ax1.set_ylim(-60.,60.)
        ax1.set_yticks(np.arange(-60.,61.,30.))
        ax1.set_xticks(np.arange(0.,361.,90.))
        ax1.set_title(title)
        land_mask = '/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc'
        land = xr.open_dataset(land_mask)
        land.lsm[0,:,:].plot.contour(ax=ax1, x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        
        ax1.set_ylabel('Latitude')
        ax1.set_xlabel('Longitude')
        
        plt.subplots_adjust(left=0.1, right=0.97, top=0.93, bottom=0.05, hspace=0.25, wspace=0.2)
        cb1=fig.colorbar(f1, ax=ax1, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.15, aspect=60, shrink=0.5)
        
        windtypestr=''; msestr=''; vidstr=''
        if windtype != 'full':
            windtypestr = '_' + windtype
        if mse:
            msestr = '_mse'
        if video:
            vidstr='video/'
        
        plot_dir = '/scratch/rg419/plots/zonal_asym_runs/gill_development/jra/' + vidstr + windtype + msestr + '/' 
        mkdir = sh.mkdir.bake('-p')
        mkdir(plot_dir)
        
        if video:
            plt.savefig(plot_dir + 'wind_and_slp_zanom_' + str(int(data_u.pentad[i])) + windtypestr + msestr + '.png', format='png')
        else:
            plt.savefig(plot_dir + 'wind_and_slp_zanom_' + str(int(data_u.pentad[i])) + windtypestr + msestr + '.pdf', format='pdf')
        plt.close()
        

import subprocess    
def make_video(filepattern, output, startno=20):
    
    command = 'ffmpeg  -framerate 5 -y -start_number '+str(startno)+' -i ' + filepattern + ' -vframes 45 -c:v libx264 -r 6 -pix_fmt yuv420p -vf scale=3200:-2 ' + output    
    subprocess.call([command], shell=True)

if __name__ == "__main__":
    
    plot_gill_dev(mse=True)
    #plot_gill_dev(windtype='none', video=True)
    
    #plot_gill_dev(windtype='rot')
    #plot_gill_dev(windtype='div')
    
    
    #make_video('/scratch/rg419/plots/zonal_asym_runs/gill_development/jra/video/none/wind_and_slp_zanom_%02d_none.png', 
    #                  '/scratch/rg419/plots/zonal_asym_runs/gill_development/jra/video/none/precip_and_slp_anom_sh.mp4', startno=56)
    
                            
                            