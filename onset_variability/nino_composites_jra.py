'''04/04/2019 Make plots of JRA-55 data fields for: March nino3.4 > or < 0, and nino3.4 < -0.5, > 0.5, and between +-0.5
Also test Jan and Feb'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
from mpl_toolkits.mplot3d import Axes3D
import statsmodels.api as sm
from data_handling_updates import gradients as gr, model_constants as mc
import pandas as pd

plot_dir = '/scratch/rg419/plots/onset_variability_new/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

# Load in nino_34 index
nino_34 = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/nino_34.nc')
nino_34_march = nino_34.nino_34.sel(time=pd.to_datetime(['%04d-03-15' % y for y in range(1958,2017)]))
nino_34_jan = nino_34.nino_34.sel(time=pd.to_datetime(['%04d-01-15' % y for y in range(1958,2017)]))
nino_34_feb = nino_34.nino_34.sel(time=pd.to_datetime(['%04d-02-15' % y for y in range(1958,2017)]))

years = np.arange(1958,2017)

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

nino_34_march = get_phase(nino_34_march)
nino_34_feb = get_phase(nino_34_feb)
nino_34_jan = get_phase(nino_34_jan)

# Plot nino_34 
month = ['March', 'Feb', 'Jan']
i=0
for nino in [nino_34_march, nino_34_feb, nino_34_jan]:
    plt.plot(nino.time,nino,'o-', color='C0')
    plt.xlabel('Time')
    plt.ylabel(month[i] + ' nino3.4 index')
    plt.yticks(np.arange(-2.,2.1,0.5))
    plt.grid(True,linestyle=':')
    plt.savefig(plot_dir + 'nino_34_test_'+ month[i] + '.pdf', format='pdf')
    plt.close()
    i=i+1


# Load in JRA data
data_u = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/ucomp_monthly/atmos_monthly_together.nc')
data_u = data_u['var33'].load().loc['1958-01':'2016-12']

# v has different time coord to u, presumably due to how Stephen has downloaded/averaged. I think the two are equivalent, so just substitute the time dimension into v
data_v_temp = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/vcomp_monthly/atmos_monthly_together.nc')
data_v = xr.DataArray(data_v_temp.var34.values, coords={'time': data_u.time, 'lev': data_v_temp.lev,
                                                  'lat': data_v_temp.lat, 'lon': data_v_temp.lon}, dims=['time','lev','lat','lon'])

dvdx = gr.ddx(data_v, a=6371.0e3) 
dudy = gr.ddy(data_u, a=6371.0e3)
data_vo = (dvdx - dudy)*86400.

data_t = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/temp_monthly/atmos_monthly_together.nc')
data_q = xr.open_dataset('/disca/share/rg419/JRA_55/sphum_monthly/atmos_monthly_together.nc')
data_z = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/height_monthly/atmos_monthly_together.nc')

data_mse = (mc.cp_air*data_t.var11 + mc.L*data_q.var51 + 9.81*data_z.var7)/1000.

land_mask = '/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc'
land = xr.open_dataset(land_mask)


def plot_winter_climate(data, title, levels_clim=np.arange(-50.,51.,5.), levels=np.arange(-5.,5.1,0.5), local=False, lev=20000.):
    
    # Add a coordinate with enso phase    
    data.coords['enso_phase'] = (('time'), np.repeat(nino_34_jan.enso_phase.values,12)) 
    data = data.sel(lev=lev)
    
    data_clim = data.groupby('time.month').mean('time').sel(month=[1,2,3]).mean('month')
    data_en = (data.where(data['enso_phase']=='El Nino', drop=True).groupby('time.month').mean('time') - data.groupby('time.month').mean('time')).sel(month=[1,2,3]).mean('month')
    data_ne= (data.where(data['enso_phase']=='Neutral', drop=True).groupby('time.month').mean('time') - data.groupby('time.month').mean('time')).sel(month=[1,2,3]).mean('month')
    data_ln = (data.where(data['enso_phase']=='La Nina', drop=True).groupby('time.month').mean('time') - data.groupby('time.month').mean('time')).sel(month=[1,2,3]).mean('month')
    
    # Start figure with 4 subplots
    rcParams['figure.figsize'] = 15, 4.5
    rcParams['font.size'] = 14
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey='row')
    axes = [ax1, ax2, ax3, ax4]
    title_list = ['Climatology', 'La Nina', 'Neutral', 'El Nino']
    
    f1 = data_clim.plot.contourf(ax=ax1, x='lon',y='lat', add_labels=False, add_colorbar=False, extend='both', levels=levels_clim)
    f2 = data_ln.plot.contourf(ax=ax2, x='lon',y='lat', add_labels=False, add_colorbar=False, extend='both', levels=levels)
    f3 = data_ne.plot.contourf(ax=ax3, x='lon',y='lat', add_labels=False, add_colorbar=False, extend='both', levels=levels)
    f4 = data_en.plot.contourf(ax=ax4, x='lon',y='lat', add_labels=False, add_colorbar=False, extend='both', levels=levels)
    
    i=0
    for ax in axes:
        land.lsm[0,:,:].plot.contour(ax=ax, x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        if local:
            ax.set_xlim(70.,150.)
            ax.set_ylim(-10.,50.)
        else:
            ax.set_xticks(np.arange(0.,361.,90.))
            ax.set_yticks(np.arange(-60.,61.,30.))
            ax.set_ylim(-60.,60.)
        ax.set_title(title_list[i])
        ax.grid(True,linestyle=':')
        ax.set_xlabel('Longitude')
        i=i+1
    
    ax1.set_ylabel('Latitude')
    
    plt.subplots_adjust(left=0.06, right=0.97, top=0.92, bottom=0.1, hspace=0.3, wspace=0.2)
    
    
    if title=='mse':
        cb1=fig.colorbar(f1, ax=ax1, use_gridspec=True, orientation = 'horizontal',fraction=0.07, pad=0.2, aspect=60, shrink=1., ticks=np.arange(250.,340.,20.))
    else:
        cb1=fig.colorbar(f1, ax=ax1, use_gridspec=True, orientation = 'horizontal',fraction=0.07, pad=0.2, aspect=30, shrink=1.)
    cb2=fig.colorbar(f2, ax=axes[1:], use_gridspec=True, orientation = 'horizontal',fraction=0.07, pad=0.2, aspect=60, shrink=0.75)
    #cb1.set_label(var)
    
    if local:
        plt.savefig('/scratch/rg419/plots/onset_variability_new/enso/JFM_enso_jan_' + title + '_local.pdf', format='pdf')
    else:
        plt.savefig('/scratch/rg419/plots/onset_variability_new/enso/JFM_enso_jan_' + title + '.pdf', format='pdf')
    plt.close()
    

plot_winter_climate(data_vo, 'vo_jra', levels_clim=np.arange(-3.,4.,1.), levels=np.arange(-0.6,0.7,0.2))
plot_winter_climate(data_vo, 'vo_jra', levels_clim=np.arange(-3.,4.,1.), levels=np.arange(-0.6,0.7,0.2), local=True)
plot_winter_climate(data_u, 'ucomp')
plot_winter_climate(data_u, 'ucomp', local=True)
plot_winter_climate(data_mse, 'mse', local=True, levels_clim=np.arange(250.,331.,10.), levels=np.arange(-2.5,2.7,0.25), lev=85000.)
plot_winter_climate(data_mse, 'mse', levels_clim=np.arange(250.,331.,10.), levels=np.arange(-2.5,2.7,0.25), lev=85000.)

data_u.close()
data_v.close()
data_vo.close()



