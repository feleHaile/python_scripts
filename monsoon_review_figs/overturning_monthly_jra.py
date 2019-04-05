''' 
25/01/2019 Plot local mass streamfunction from climatology month by month
'''
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from hadley_cell import mass_streamfunction
from data_handling_updates import model_constants as mc
from windspharm.xarray import VectorWind


def monthly_overturning_jra(lonin=[-1.,361.], region='All', sanity_check=False):   

    plot_dir = '/scratch/rg419/plots/monsoon_review_figs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    # Load in data
    data_u = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/ucomp_monthly/atmos_monthly_together.nc')
    data_v = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/vcomp_monthly/atmos_monthly_together.nc')
    
    # Make climatologies
    u_clim = data_u.groupby('time.month').mean('time')
    v_clim = data_v.groupby('time.month').mean('time')
    print('means taken')
    
    def get_lons(lonin, data):
        if lonin[1]>lonin[0]:
            lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
        else:
            lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
        return lons
    
    lons = get_lons(lonin,data_u)
    
    levs = np.arange(50.,1000.,50.)
    
    ds_clim = xr.Dataset({'ucomp': (['month', 'pfull', 'lat', 'lon'], u_clim.var33.values), 'vcomp': (['month', 'pfull', 'lat', 'lon'], v_clim.var34.values)},
                     coords={'month': ('month', v_clim.month),
                               'pfull': ('pfull', v_clim.lev/100.),
                               'lat': ('lat', v_clim.lat),
                               'lon': ('lon', v_clim.lon)})
    
    ds_clim = ds_clim.sel(pfull=levs)
                        
    psi = mass_streamfunction(ds_clim, lons=lons, dp_in=50.)
    psi /= 1.e9
    
    # Set figure parameters
    rcParams['figure.figsize'] = 14, 9
    rcParams['font.size'] = 14
    
    # Start figure with 12 subplots
    # Start figure with 12 subplots
    fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4)
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'June', 'Jul', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec']
        
    def plot_psi_u(ax, psi, u, month):
        f1 = u.sel(lon=lons, month=month).mean('lon').plot.contourf(ax=ax, x='lat', y='pfull', yincrease=False, levels=np.arange(-50.,50.1,5.), extend='both', add_labels=False, add_colorbar=False)
        
        psi.sel(month=month).plot.contour(ax=ax, x='lat', y='pfull', yincrease=False, levels=np.arange(0.,601,50.), colors='k', add_labels=False)
        psi.sel(month=month).plot.contour(ax=ax, x='lat', y='pfull', yincrease=False, levels=np.arange(-600.,0.,50.), colors='k', linestyles='dashed', add_labels=False)
        
        m = mc.omega * mc.a**2. * np.cos(u.lat*np.pi/180.)**2. + u.sel(lon=lons, month=month).mean('lon') * mc.a * np.cos(u.lat*np.pi/180.)
        m_levs = mc.omega * mc.a**2. * np.cos(np.arange(-60.,1.,5.)*np.pi/180.)**2.    
        m.plot.contour(ax=ax, x='lat', y='pfull', yincrease=False, levels=m_levs, colors='0.7', add_labels=False)
        
        ax.set_xlim(-35,35)
        ax.set_xticks(np.arange(-30,31,15))
        ax.grid(True,linestyle=':')
        ax.set_title(months[i-1])
        
        return f1
    
    for i in range(1,13):
        f1 = plot_psi_u(axes[i-1], psi, ds_clim.ucomp, i)
    
    for i in range(8,12):
        axes[i].set_xlabel('Latitude')
    for i in range(0,12,4):
        axes[i].set_ylabel('Pressure, hPa')
        
    plt.subplots_adjust(left=0.1, right=0.97, top=0.95, bottom=0.1, hspace=0.3, wspace=0.3)
    
    cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.1, aspect=30, shrink=0.5)
    cb1.set_label('Zonal wind speed, m/s')
    
    if lonin == [-1.,361.]:
        plt.savefig(plot_dir + 'monthly_psi_u_jra.pdf', format='pdf')
    else:
        figname = plot_dir + 'monthly_psi_u_jra_' + region + '.pdf'
        plt.savefig(figname, format='pdf')

        
    plt.close()
    
    
#monthly_overturning_jra()
monthly_overturning_jra([60.,80.], 'India')
monthly_overturning_jra([90.,110.], 'BoB')
monthly_overturning_jra([110.,140.], 'East Asia')
monthly_overturning_jra([140.,160.], 'Western North Pacific')
monthly_overturning_jra([240.,280.], 'North America')
monthly_overturning_jra([355.,5.], 'West Africa')
monthly_overturning_jra([25.,45.], 'South Africa')
monthly_overturning_jra([300.,320.], 'South America')
