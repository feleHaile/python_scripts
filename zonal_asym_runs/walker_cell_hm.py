''' 
16/08/2018 Plot hm of local walker circulation
'''
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from hadley_cell import walker_cell
from data_handling_updates import model_constants as mc
from windspharm.xarray import VectorWind


def walker_hm(run, regions=[[0,10], [10,20], [20,30], [30,40]]):
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
        
    plot_dir = '/scratch/rg419/plots/overturning_monthly/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    # Create a VectorWind instance to handle the computation
    w = VectorWind(data.ucomp.sel(pfull=np.arange(50.,950.,50.)), data.vcomp.sel(pfull=np.arange(50.,950.,50.)))
    uchi, vchi, upsi, vpsi = w.helmholtz()
    
    # Set figure parameters
    rcParams['figure.figsize'] = 10, 7
    rcParams['font.size'] = 14
    # Start figure with 4 subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
    axes = [ax1,ax2,ax3,ax4]
    
    i=0
    for ax in axes:
        psi_w = walker_cell(uchi, latin=regions[i], dp_in=-50.)
        psi_w /= 1.e9
        i=i+1
        #f1=psi_w.sel(pfull=500).plot.contourf(ax=ax, x='lon', y='xofyear', add_labels=False, add_colorbar=False, levels=np.arange(-500.,501.,100.), extend='both')
        f1=psi_w.sel(pfull=500).plot.contourf(ax=ax, x='lon', y='xofyear', add_labels=False, add_colorbar=False, levels=np.arange(-200.,201.,20.), extend='both')
        
    ax1.set_title('0-10 N')
    ax2.set_title('10-20 N')
    ax3.set_title('20-30 N')
    ax4.set_title('30-40 N')
    
    for ax in [ax1,ax2,ax3,ax4]:
        ax.grid(True,linestyle=':')
        ax.set_yticks(np.arange(0.,73., 18.))
        ax.set_xticks([0.,90.,180.,270.,360.])
    
    ax3.set_xlabel('Longitude')
    ax4.set_xlabel('Longitude')
    ax1.set_ylabel('Pentad')
    ax3.set_ylabel('Pentad')
    
    plt.subplots_adjust(left=0.1, right=0.97, top=0.95, bottom=0.1, hspace=0.3, wspace=0.3)
    
    cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.1, aspect=30, shrink=0.5)
    cb1.set_label('Zonal overturning streamfunction')
    
    # Save as a pdf
    plt.savefig(plot_dir + 'walker_cell_hm_' + run + '.pdf', format='pdf')
    plt.close()
    
    data.close()

walker_hm('half_shallow')
#walker_hm('half_nh_shallow')
#walker_hm('half_10_shallow')
#walker_hm('half_30_shallow')
#walker_hm('half_shallow_5')
#walker_hm('half_shallow_10')
#walker_hm('3q_shallow')
#walker_hm('q_shallow')

