''' 
07/09/2018 Plot local mass vertical mass flux associated with divergent flow in Hadley and Walker cells in NH winter and summer
'''
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from hadley_cell import mass_streamfunction
from data_handling_updates import model_constants as mc, gradients as gr
from windspharm.xarray import VectorWind


def h_w_mass_flux_monthly(lev=50000., dp=5000.):
    
    # Load in data
    data_u = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/ucomp_monthly/atmos_monthly_together.nc')
    data_v = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/vcomp_monthly/atmos_monthly_together.nc')
    land_mask='/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc'
    land = xr.open_dataset(land_mask)
    
    data_u = data_u.sel(time=slice('1979','2016'))
    data_v = data_v.sel(time=slice('1979','2016'))
    
    # Make climatologies
    u_clim = data_u.groupby('time.month').mean('time')
    v_clim = data_v.groupby('time.month').mean('time')
    print('means taken')
    
    plot_dir = '/scratch/rg419/plots/monsoon_review_figs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    # Create a VectorWind instance to handle the computation
    w = VectorWind(u_clim.var33, v_clim.var34)
    # Compute variables
    streamfun, vel_pot = w.sfvp()
    uchi, vchi, upsi, vpsi = w.helmholtz()
    coslat = np.cos(data_u.lat * np.pi/180)
    
    # Evaluate mass fluxes for the zonal and meridional components (Walker and Hadley) following Schwendike et al. 2014
    mass_flux_zon = (gr.ddx(uchi.sel(lev=np.arange(5000.,100000., 5000.)))).cumsum('lev') * dp * coslat/ mc.grav
    mass_flux_merid = (gr.ddy(vchi.sel(lev=np.arange(5000.,100000., 5000.)))).cumsum('lev') * dp * coslat/ mc.grav
    
    # Set figure parameters
    rcParams['figure.figsize'] = 15, 9
    rcParams['font.size'] = 16
    
    # Start figure with 4 subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
    axes = [ax1, ax2, ax3, ax4]

    f1 = mass_flux_merid.sel(lev=lev, month=[5,6,7,8,9]).mean('month').plot.contourf(
              ax=ax1, x='lon', y='lat', add_labels=False, add_colorbar=False, levels=np.arange(-0.006,0.007,0.001), extend='both')
    f1 = mass_flux_zon.sel(lev=lev, month=[5,6,7,8,9]).mean('month').plot.contourf(
              ax=ax2, x='lon', y='lat', add_labels=False, add_colorbar=False, levels=np.arange(-0.006,0.007,0.001), extend='both')
    f1 = mass_flux_merid.sel(lev=lev, month=[11,12,1,2,3]).mean('month').plot.contourf(
              ax=ax3, x='lon', y='lat', add_labels=False, add_colorbar=False, levels=np.arange(-0.006,0.007,0.001), extend='both')
    f1 = mass_flux_zon.sel(lev=lev, month=[11,12,1,2,3]).mean('month').plot.contourf(
              ax=ax4, x='lon', y='lat', add_labels=False, add_colorbar=False, levels=np.arange(-0.006,0.007,0.001), extend='both')
    
    i=0
    for ax in axes:
        land.lsm[0,:,:].plot.contour(ax=ax, x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', linewidths=2.)
        ax.fill_between([0,360], [-30, -30], [-10, -10], alpha=0.1, color='k')
        ax.fill_between([0,360], [10, 10], [30, 30], alpha=0.1, color='k')
        ax.set_ylim(-60,60)
        ax.set_xticks(np.arange(0,361,90))
        ax.set_yticks(np.arange(-60,61,30))
        ax.grid(True,linestyle=':')
        ax.set_xlim(0,360)
    
    ax1.text(-40, 60, 'a)')
    ax2.text(-30, 60, 'b)')
    ax3.text(-40, 60, 'c)')
    ax4.text(-30, 60, 'd)')
    plt.subplots_adjust(left=0.05, right=0.97, top=0.95, bottom=0.05, hspace=0.1, wspace=0.1)
    
    cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.1, aspect=30, shrink=0.5)
    cb1.set_label('Vertical mass flux, kgm$^{-2}$s$^{-1}$')
    
    figname = plot_dir + 'hadley_walker_jra.pdf'
    plt.savefig(figname, format='pdf')
    plt.close()
    

h_w_mass_flux_monthly()
