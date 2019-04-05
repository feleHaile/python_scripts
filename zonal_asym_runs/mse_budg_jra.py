'''15/3/2019 Plot single level MSE budget for reanalysis
'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from data_handling_updates import gradients as gr
import matplotlib.patches as patches

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

def mse_budg_jra():
    
    cp_air = 1004.6
    L = 2.507e6
    grav = 9.80665
    
    plot_dir = '/scratch/rg419/plots/zonal_asym_runs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data = xr.open_dataset('/disca/share/rg419/JRA_55/heating_terms.nc')
#    data_w = xr.open_dataset('/disca/share/rg419/jra_omega_pentad_clim_alllevs.nc')
    
    mse = (cp_air * data.temp + L * data.sphum + grav * data.height)
    dse = (cp_air * data.temp + grav * data.height)
    qse = (L * data.sphum)

    q_diab = L * (data.cnvmr + data.lrgmr + data.vdfmr)/86400.
    t_diab = cp_air * (data.cnvhr + data.lrghr + data.vdfhr + data.lwhr + data.swhr)/86400.
    q_total = t_diab + q_diab
    
    u_mse = (cp_air * data.ucomp_temp + L * data.sphum_u + grav * data.ucomp_height)
    v_mse = (cp_air * data.vcomp_temp + L * data.sphum_v + grav * data.vcomp_height)
    w_mse = (cp_air * data.omega_temp + L * data.sphum_w + grav * data.omega_height)
    u_dse = (cp_air * data.ucomp_temp + grav * data.ucomp_height)
    v_dse = (cp_air * data.vcomp_temp + grav * data.vcomp_height)
    w_dse = (cp_air * data.omega_temp + grav * data.omega_height)
    u_qse = (L * data.sphum_u)
    v_qse = (L * data.sphum_v)
    w_qse = (L * data.sphum_w)
    
    def get_budget_terms(energy, uenergy, venergy, wenergy, diab):
        u_e_eddy = uenergy - energy * data.ucomp
        v_e_eddy = venergy - energy * data.vcomp
        w_e_eddy = wenergy - energy * data.omega
        eddy_conv = -1.*(gr.ddx(u_e_eddy) + gr.ddy(v_e_eddy) + gr.ddp(w_e_eddy))
        horiz_adv = -1.*(data.ucomp * gr.ddx(energy) + data.vcomp * gr.ddy(energy, vector=False))
        denergydt = gr.ddt(energy)*20.
        total = (diab + horiz_adv + vert_adv + eddy_conv)*10.
        return diab, horiz_adv, vert_adv, eddy_conv, denergydt, total
    
    mse_terms = get_budget_terms(mse, u_mse, v_mse, w_mse, q_total)
    dse_terms = get_budget_terms(dse, u_dse, v_dse, w_dse, t_diab)
    qse_terms = get_budget_terms(qse, u_qse, v_qse, w_qse, q_diab)
    
    
    
    for pentad in [32,38,44,50,56]:
        # Start figure with 4 subplots
        rcParams['figure.figsize'] = 17, 7
        rcParams['font.size'] = 14
        fig, ((ax1, ax2, ax3, ax4, ax5, ax6), (ax7, ax8, ax9, ax10, ax11, ax12), (ax13, ax14, ax15, ax16, ax17, ax18)) = plt.subplots(3, 6, sharex='col', sharey='row')
        axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12, ax13, ax14, ax15, ax16, ax17, ax18]
        
        i=0
        for energy in [mse_terms, dse_terms, qse_terms]:
            f1 = energy[0].sel(xofyear=pentad).plot.contourf(ax=axes[i*6],x='lon',y='lat', add_labels=False, extend='both', levels=np.arange(-0.1, 0.11, 0.01), add_colorbar=False)
            energy[1].sel(xofyear=pentad).plot.contourf(ax=axes[i*6+1],x='lon',y='lat', add_labels=False,  extend='both', levels=np.arange(-0.1, 0.11, 0.01), add_colorbar=False)
            energy[2].sel(xofyear=pentad).plot.contourf(ax=axes[i*6+2],x='lon',y='lat', add_labels=False, extend='both', levels=np.arange(-0.1, 0.11, 0.01), add_colorbar=False)
            energy[3].sel(xofyear=pentad).plot.contourf(ax=axes[i*6+3],x='lon',y='lat', add_labels=False, extend='both', levels=np.arange(-0.1, 0.11, 0.01), add_colorbar=False)
            energy[4].sel(xofyear=pentad).plot.contourf(ax=axes[i*6+4],x='lon',y='lat', add_labels=False, extend='both', levels=np.arange(-0.1, 0.11, 0.01), add_colorbar=False)
            energy[5].sel(xofyear=pentad).plot.contourf(ax=axes[i*6+5],x='lon',y='lat', add_labels=False, extend='both', levels=np.arange(-0.1, 0.11, 0.01), add_colorbar=False)
            i=i+1
        
        for ax in axes:
            ax.set_ylim(-15.,45.)
            ax.set_xlim(90.,225.)
            ax.set_yticks(np.arange(-15.,45.,15.))
            ax.set_xticks(np.arange(90.,226.,45.))
            ax.grid(True,linestyle=':')
        
        land_mask = '/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc'
        land = xr.open_dataset(land_mask)
        land.lsm[0,:,:].plot.contour(ax=ax, x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        (land.z[0,:,:]/9.81).plot.contour(ax=ax, x='longitude', y='latitude', levels=np.arange(2000.,3001.,1000.), add_labels=False, colors='k')
    
        plt.subplots_adjust(left=0.1, right=0.97, top=0.93, bottom=0.05, hspace=0.25, wspace=0.2)
        cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.05, aspect=60, shrink=1.)
        
        plt.savefig(plot_dir + 'mse_budg_850_JRA55_' + str(pentad) + '.pdf', format='pdf')
        plt.close()
        

if __name__ == "__main__":
    
    mse_budg_jra()
                          
                            