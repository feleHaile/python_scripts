'''12/3/2019 Plot vertically integrated MSE budget
'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from data_handling_updates import gradients as gr, model_constants as mc
import matplotlib.patches as patches


def mse_budg_vert_int(run, land_mask=None):
    
    plot_dir = '/scratch/rg419/plots/zonal_asym_runs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    data['rflux_surf'] = data.t_surf ** 4. * mc.stefan - data.flux_sw - data.flux_lw
    data['rflux_toa'] = data.toa_sw - data.olr
    
    mse = (mc.cp_air * data.temp + mc.L * data.sphum + mc.grav * data.height)/1000.
    e = ((mc.cp_air - mc.rdgas) * data.temp + mc.L * data.sphum + mc.grav * data.height)/1000.
    u_mse = (mc.cp_air * data.ucomp_temp + mc.L * data.sphum_u + mc.grav * data.ucomp_height)/1000.
    v_mse = (mc.cp_air * data.vcomp_temp + mc.L * data.sphum_v + mc.grav * data.vcomp_height)/1000.
    w_mse = (mc.cp_air * data.omega_temp + mc.L * data.sphum_w + mc.grav * data.omega_height)/1000.
    
    u_mse_eddy = u_mse - mse * data.ucomp
    v_mse_eddy = v_mse - mse * data.vcomp
    w_mse_eddy = w_mse - mse * data.omega
    
    u_dmsedx = data.ucomp * gr.ddx(mse)
    v_dmsedy = data.vcomp * gr.ddy(mse, vector=False)
    w_dmsedp = data.omega * gr.ddp(mse)
    dumsedx_eddy = gr.ddx(u_mse_eddy)
    dvmsedy_eddy = gr.ddy(v_mse_eddy)
    dwmsedp_eddy = gr.ddp(w_mse_eddy)
    dedt = gr.ddt(e)*5.
    
    Fnet = data.flux_lhe + data.flux_t + data.rflux_surf + data.rflux_toa
    
    def column_int(var_in):
        var_int = mc.cp_air * var_in.sum('pfull')*5000./mc.grav
        return var_int
    
    u_dmsedx = column_int(u_dmsedx)
    v_dmsedy = column_int(v_dmsedy)
    w_dmsedp = -1.*column_int(w_dmsedp)
    dumsedx_eddy = column_int(dumsedx_eddy)
    dvmsedy_eddy = column_int(dvmsedy_eddy)
    dwmsedp_eddy = column_int(dwmsedp_eddy)
    dedt = column_int(dedt)
    
    horiz_adv = -1.*(u_dmsedx + v_dmsedy)
    eddy_conv = -1.*(dumsedx_eddy + dvmsedy_eddy + dwmsedp_eddy)
    
    total = (Fnet + horiz_adv + eddy_conv + w_dmsedp)*5.
    
    for pentad in [32,38,44,50,56]:
        # Start figure with 4 subplots
        rcParams['figure.figsize'] = 14, 5
        rcParams['font.size'] = 11
        fig, ((ax1, ax2, ax3, ax4, ax5), (ax6, ax7, ax8, ax9, ax10)) = plt.subplots(2, 5, sharex='col', sharey='row')
        axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10]
        
        f1 = Fnet.sel(xofyear=pentad).plot.contourf(ax=ax1,x='lon',y='lat', add_labels=False, extend='both', levels=np.arange(-200.,201.,20.), add_colorbar=False)
        horiz_adv.sel(xofyear=pentad).plot.contourf(ax=ax2,x='lon',y='lat', add_labels=False,  extend='both', levels=np.arange(-200.,201.,20.), add_colorbar=False)
        w_dmsedp.sel(xofyear=pentad).plot.contourf(ax=ax3,x='lon',y='lat', add_labels=False, extend='both', levels=np.arange(-200.,201.,20.), add_colorbar=False)
        eddy_conv.sel(xofyear=pentad).plot.contourf(ax=ax4,x='lon',y='lat', add_labels=False, extend='both', levels=np.arange(-200.,201.,20.), add_colorbar=False)
        dedt.sel(xofyear=pentad).plot.contourf(ax=ax5,x='lon',y='lat', add_labels=False, extend='both', levels=np.arange(-200.,201.,20.), add_colorbar=False)
        
        data.flux_lhe.sel(xofyear=pentad).plot.contourf(ax=ax6,x='lon',y='lat', add_labels=False, extend='both', levels=np.arange(-200.,201.,20.), add_colorbar=False)
        data.flux_t.sel(xofyear=pentad).plot.contourf(ax=ax7,x='lon',y='lat', add_labels=False, extend='both', levels=np.arange(-200.,201.,20.), add_colorbar=False)
        data.rflux_surf.sel(xofyear=pentad).plot.contourf(ax=ax8,x='lon',y='lat', add_labels=False, extend='both', levels=np.arange(-200.,201.,20.), add_colorbar=False)
        data.rflux_toa.sel(xofyear=pentad).plot.contourf(ax=ax9,x='lon',y='lat', add_labels=False, extend='both', levels=np.arange(-200.,201.,20.), add_colorbar=False)
        total.sel(xofyear=pentad).plot.contourf(ax=ax10,x='lon',y='lat', add_labels=False, extend='both', levels=np.arange(-200.,201.,20.), add_colorbar=False)
        
        
        for ax in axes:
            ax.set_ylim(-15.,45.)
            ax.set_xlim(90.,225.)
            ax.set_yticks(np.arange(-15.,45.,15.))
            ax.set_xticks(np.arange(90.,226.,45.))
            ax.grid(True,linestyle=':')
        
        
        if not land_mask==None:
            land = xr.open_dataset(land_mask)
            for ax in axes:
                land.land_mask.plot.contour(x='lon', y='lat', ax=ax, levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        
        plt.subplots_adjust(left=0.05, right=0.97, top=0.97, bottom=0.03, hspace=0.25, wspace=0.2)
        cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.05, aspect=60, shrink=1.)
        
        plt.savefig(plot_dir + 'mse_budg_' + run + '_' + str(pentad) + '.pdf', format='pdf')
        plt.close()
        

if __name__ == "__main__":
    
    mse_budg_vert_int('half_land', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')
    mse_budg_vert_int('half_nh_land_tibet_ol8', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_nh_land_tibet.nc')
                          
                            