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


def plot_vort_dev(video=False, threed=True):
    
    data_w = xr.open_dataset('/disca/share/rg419/jra_omega_pentad_clim.nc')
    data_u = xr.open_dataset('/disca/share/rg419/jra_ucomp_pentad_clim.nc')
    data_v_temp = xr.open_dataset('/disca/share/rg419/jra_vcomp_pentad_clim.nc')
    # v has different time coord to u, presumably due to how Stephen has downloaded/averaged. I think the two are equivalent, so just substitute the time dimension into v
    data_v = xr.DataArray(data_v_temp.var34.values, coords={'pentad': data_u.pentad, 'lat': data_u.lat, 'lon': data_u.lon}, dims=('pentad','lat','lon'))    
    print('files opened')
    
    data_u = data_u.var33
    data_w = data_w.var39
    
    zon_adv = data_u.sel(lev=20000.) * gr.ddx(data_u.sel(lev=20000.))
    merid_adv = data_v * gr.ddy(data_u.sel(lev=20000.))
    vert_adv = data_w * (gr.ddp(data_u, pname='lev')).sel(lev=20000.)*100.
    
    sinphi = np.sin(data_u.lat * np.pi/180.)
    f = 2.* mc.omega*sinphi
    if threed:
        rossby =  (zon_adv + merid_adv + vert_adv)/(f*data_v)
    else:
        rossby = merid_adv/(f*data_v)
    # Start figure with 1 subplots
    rcParams['figure.figsize'] = 10, 5
    rcParams['font.size'] = 14

    for i in range(72):
        fig, ax1 = plt.subplots()
        title = 'Pentad ' + str(int(data_u.pentad[i]))
        
        f1 = rossby.sel(pentad=i+1).plot.contourf(x='lon', y='lat', ax=ax1, add_labels=False, add_colorbar=False, extend='both', zorder=1, levels = np.arange(0.,1.1,0.1))
        
#        data_p.precipitation.sel(pentad=i+1).plot.contour(x='lon', y='lat', ax=ax1, add_labels=False, extend='both', zorder=1, levels=np.arange(2.,21.,2.), colors='k')        
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
    
        vidstr=''
        if video:
            vidstr='video/'
        
        if threed:
            plot_dir = '/scratch/rg419/plots/zonal_asym_runs/gill_development/jra/' + vidstr + '/rossby3d/' 
        else:
            plot_dir = '/scratch/rg419/plots/zonal_asym_runs/gill_development/jra/' + vidstr + '/rossby/' 
        mkdir = sh.mkdir.bake('-p')
        mkdir(plot_dir)
        
        if video:
            plt.savefig(plot_dir + 'rossby_and_precip_' + str(int(data_u.pentad[i])) + '.png', format='png')
        else:
            plt.savefig(plot_dir + 'rossby_and_precip_' + str(int(data_u.pentad[i])) + '.pdf', format='pdf')
        plt.close()


import subprocess
def make_video(filepattern, output):
    command = 'ffmpeg  -framerate 5 -y -start_number 30 -i ' + filepattern + ' -vframes 45 -c:v libx264 -r 6 -pix_fmt yuv420p -vf scale=3200:-2 ' + output    
    subprocess.call([command], shell=True)
    

if __name__ == "__main__":

    #plot_vort_dev(threed=True)
    plot_vort_dev(threed=False)
   # plot_vort_dev('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc', video=True, threed=True)

    #make_video('/scratch/rg419/plots/zonal_asym_runs/gill_development/half_shallow/video/rossby3d/rossby_and_precip_%02d.png', 
    #                  '/scratch/rg419/plots/zonal_asym_runs/gill_development/half_shallow/video/rossby3d/rossby_and_precip.mp4')

    
                         