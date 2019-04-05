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


def summer_overturning_jra(lonin=[-1.,361.], region='All', sanity_check=False):   

    plot_dir = '/scratch/rg419/plots/monsoon_review_figs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    # Load in data
    data_u = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/ucomp_monthly/atmos_monthly_together.nc')
    data_v = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/vcomp_monthly/atmos_monthly_together.nc')
    
    # Make climatologies
    u_clim = data_u.groupby('time.month').mean('time')
    v_clim = data_v.groupby('time.month').mean('time')
    
    summer_a = [11,12,1,2,3]
    summer_b = [5,6,7,8,9]
    
    # Take austral and boreal summer means
    u_summer_a = u_clim.sel(month=summer_a).mean('month') 
    v_summer_a = v_clim.sel(month=summer_a).mean('month')
    
    u_summer_b = u_clim.sel(month=summer_b).mean('month') 
    v_summer_b = v_clim.sel(month=summer_b).mean('month')
    print('means taken')
    
    def get_lons(lonin, data):
        if lonin[1]>lonin[0]:
            lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
        else:
            lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
        return lons
    
    lons = get_lons(lonin,data_u)
    
    levs = np.arange(50.,1000.,50.)
    
    ds_summer_a = xr.Dataset({'ucomp': (['pfull', 'lat', 'lon'], u_summer_a.var33.values), 'vcomp': (['pfull', 'lat', 'lon'], v_summer_a.var34.values)},
                     coords={'pfull': ('pfull', v_summer_a.lev/100.),
                               'lat': ('lat', v_summer_a.lat),
                               'lon': ('lon', v_summer_a.lon)})
    
    ds_summer_b = xr.Dataset({'ucomp': (['pfull', 'lat', 'lon'], u_summer_b.var33.values), 'vcomp': (['pfull', 'lat', 'lon'], v_summer_b.var34.values)},
                     coords={'pfull': ('pfull', v_summer_a.lev/100.),
                               'lat': ('lat', v_summer_a.lat),
                               'lon': ('lon', v_summer_a.lon)})
    
    ds_summer_a = ds_summer_a.sel(pfull=levs)
    ds_summer_b = ds_summer_b.sel(pfull=levs)
                        
    psi_summer_a = mass_streamfunction(ds_summer_a, lons=lons, dp_in=50.)
    psi_summer_a /= 1.e9
    
    psi_summer_b = mass_streamfunction(ds_summer_b, lons=lons, dp_in=50.)
    psi_summer_b /= 1.e9
            

    # Set figure parameters
    rcParams['figure.figsize'] = 5, 7
    rcParams['font.size'] = 14
    
    # Start figure with 12 subplots
    fig, (ax1, ax2) = plt.subplots(2, 1)
    
    def plot_psi_u(ax, psi, u):
        f1 = u.sel(lon=lons).mean('lon').plot.contourf(ax=ax, x='lat', y='pfull', yincrease=False, levels=np.arange(-50.,50.1,5.), extend='both', add_labels=False, add_colorbar=False)
        
        psi.plot.contour(ax=ax, x='lat', y='pfull', yincrease=False, levels=np.arange(0.,601,50.), colors='k', add_labels=False)
        psi.plot.contour(ax=ax, x='lat', y='pfull', yincrease=False, levels=np.arange(-600.,0.,50.), colors='k', linestyles='dashed', add_labels=False)
        
        m = mc.omega * mc.a**2. * np.cos(u.lat*np.pi/180.)**2. + u.sel(lon=lons).mean('lon') * mc.a * np.cos(u.lat*np.pi/180.)
        m_levs = mc.omega * mc.a**2. * np.cos(np.arange(-60.,1.,5.)*np.pi/180.)**2.    
        m.plot.contour(ax=ax, x='lat', y='pfull', yincrease=False, levels=m_levs, colors='0.7', add_labels=False)
        
        ax.set_xlim(-35,35)
        ax.set_xticks(np.arange(-30,31,15))
        ax.grid(True,linestyle=':')
        return f1
    
    f1 = plot_psi_u(ax1, psi_summer_b, ds_summer_b.ucomp)
    plot_psi_u(ax2, psi_summer_a, ds_summer_a.ucomp)
    
    ax1.set_title(region)
    
    plt.subplots_adjust(left=0.1, right=0.97, top=0.95, bottom=0.1, hspace=0.3, wspace=0.3)
    
    cb1=fig.colorbar(f1, ax=[ax1,ax2], use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.1, aspect=30, shrink=0.5)
    cb1.set_label('Zonal wind speed, m/s')
    
    if lonin == [-1.,361.]:
        plt.savefig(plot_dir + 'summer_psi_u_jra.pdf', format='pdf')
    else:
        figname = plot_dir + 'summer_psi_u_jra_' + region + '.pdf'
        plt.savefig(figname, format='pdf')
        
    ax1.set_ylabel('Pressure, hPa')
    ax2.set_ylabel('Pressure, hPa')
    ax2.set_xlabel('Latitude')
        
    plt.close()
    
    if sanity_check:
        plt.figure(1)
        ds_summer_b.vcomp.sel(lon=lons).mean('lon').plot.contourf(x='lat', y='pfull', yincrease=False)
        plt.xlim(-35,35)
        plt.figure(2)
        ds_summer_a.vcomp.sel(lon=lons).mean('lon').plot.contourf(x='lat', y='pfull', yincrease=False)
        plt.xlim(-35,35)
        plt.show()
    
#summer_overturning_jra()
summer_overturning_jra([60.,80.], 'India')
summer_overturning_jra([90.,110.], 'BoB')
summer_overturning_jra([110.,140.], 'East Asia')
summer_overturning_jra([140.,160.], 'Western North Pacific')
summer_overturning_jra([240.,280.], 'North America')
summer_overturning_jra([355.,5.], 'West Africa')
summer_overturning_jra([25.,45.], 'South Africa')
summer_overturning_jra([300.,320.], 'South America')
