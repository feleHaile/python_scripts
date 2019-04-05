"""
Plot precip and mse hm and overturning/ascent for 10m run

"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from climatology import precip_mse_plot
from pylab import rcParams
from hadley_cell import mass_streamfunction
from data_handling_updates import make_sym

plot_dir = '/scratch/rg419/plots/paper_2_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
    
rcParams['figure.figsize'] = 5.5, 4.2
rcParams['font.size'] = 14

#fig, ax1 = plt.subplots(111)
fig = plt.figure()
ax1 = fig.add_subplot(111)    

def precip_psi_plot(run, ax):
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    psi = mass_streamfunction(data, a=6376.0e3, dp_in=50.)
    psi /= 1.e9
    
    psi = make_sym(psi, asym=True)
    data['precipitation'] = make_sym(data.precipitation)
    
    f1 = precip_mse_plot(data, ax, plot_type='precip', precip_contour=None)
    
    psi.sel(pfull=500).plot.contour(ax=ax, x='xofyear', y='lat', levels=np.arange(-500.,0.,100.), add_labels=False, colors='0.7', linewidths=2, linestyles='--')
    psi.sel(pfull=500).plot.contour(ax=ax, x='xofyear', y='lat', levels=np.arange(0.,510.,100.), add_labels=False, colors='0.7', linewidths=2)
    psi.sel(pfull=500).plot.contour(ax=ax, x='xofyear', y='lat', levels=np.arange(-1000.,1010.,1000.), add_labels=False, colors='0.5', linewidths=2)
    
    return f1

f1 = precip_psi_plot('sn_1_sst_zs', ax=ax1)
#precip_psi_plot('sine_sst_10m_zs', ax=ax2)

ax1.set_xlabel('Pentad')

#for ax in [ax1, ax2]:
ax1.set_xticks([12.,24.,36.,48.,60.,72.])
ax1.set_xlim([1,72])
ax1.fill_between([32.,35.], -60, 60,
                facecolor='k', alpha=0.3)
ax1.fill_between([39.,42.], -60, 60,
                facecolor='k', alpha=0.3)
ax1.fill_between([46.,49.], -60, 60,
                facecolor='k', alpha=0.3)

plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.07)

cb1=fig.colorbar(f1, ax=ax1, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.2, aspect=40)
cb1.set_label('Precipitation, mm/day')

plt.savefig(plot_dir+'precip_mse_hm_10m_zs.pdf', format='pdf')
plt.close()        
