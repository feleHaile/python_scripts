''' 
15/02/2019 Plot local mass streamfunction for boreal winter and summer in each monsoon region
'''
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from hadley_cell import mass_streamfunction
from data_handling_updates import model_constants as mc


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


def plot_overturning(data_u, data_v, ax, lonin=[-1.,361.], region='All', sanity_check=False):   
    
    lons = get_lons(lonin,data_u)
    levs = np.arange(50.,1000.,50.)
    
    ds_summer = xr.Dataset({'ucomp': (['pfull', 'lat', 'lon'], data_u.var33.values), 'vcomp': (['pfull', 'lat', 'lon'], data_v.var34.values)},
                     coords={'pfull': ('pfull', v_summer_a.lev/100.),
                               'lat': ('lat', v_summer_a.lat),
                               'lon': ('lon', v_summer_a.lon)})
    
    ds_summer = ds_summer.sel(pfull=levs)
                        
    psi = mass_streamfunction(ds_summer, lons=lons, dp_in=50.)#, use_v_locally=True)
    psi /= 1.e9
    
    u = ds_summer.ucomp
    f1 = u.sel(lon=lons).mean('lon').plot.contourf(ax=ax, x='lat', y='pfull', yincrease=False, levels=np.arange(-50.,50.1,5.), extend='both', add_labels=False, add_colorbar=False)
    
    psi.plot.contour(ax=ax, x='lat', y='pfull', yincrease=False, levels=np.arange(0.,601,50.), colors='k', add_labels=False)
    psi.plot.contour(ax=ax, x='lat', y='pfull', yincrease=False, levels=np.arange(-600.,0.,50.), colors='k', linestyles='dashed', add_labels=False)
    
    m = mc.omega * mc.a**2. * np.cos(u.lat*np.pi/180.)**2. + u.sel(lon=lons).mean('lon') * mc.a * np.cos(u.lat*np.pi/180.)
    m_levs = mc.omega * mc.a**2. * np.cos(np.arange(-60.,1.,5.)*np.pi/180.)**2.    
    m.plot.contour(ax=ax, x='lat', y='pfull', yincrease=False, levels=m_levs, colors='k', alpha=0.7, add_labels=False)
    
    ax.set_xlim(-35,35)
    ax.set_xticks(np.arange(-30,31,15))
    ax.grid(True,linestyle=':')
    
    ax.set_title(region)
    
    plt.subplots_adjust(left=0.05, right=0.97, top=0.95, bottom=0.05, hspace=0.3, wspace=0.2)
    
    return f1        
    


# Set figure parameters
rcParams['figure.figsize'] = 15, 6
rcParams['font.size'] = 14

# Start figure with 12 subplots
fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4, sharex='col', sharey='row')
axes=[ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
plot_overturning(u_summer_b, v_summer_b, ax1, [70.,100.], 'South Asia')
plot_overturning(u_summer_b, v_summer_b, ax2, [110.,140.], 'East Asia')
plot_overturning(u_summer_b, v_summer_b, ax3, [140.,160.], 'Western North Pacific')
plot_overturning(u_summer_a, v_summer_a, ax4, [100.,170.], 'Australia')
plot_overturning(u_summer_b, v_summer_b, ax5, [325.,45.], 'West Africa')
plot_overturning(u_summer_a, v_summer_a, ax6, [20.,90.], 'South Africa')
plot_overturning(u_summer_b, v_summer_b, ax7, [240.,280.], 'North America')
f1 = plot_overturning(u_summer_a, v_summer_a, ax8, [280.,320.], 'South America')
  
ax1.set_ylabel('Pressure, hPa')  
ax5.set_ylabel('Pressure, hPa')
for ax in [ax5, ax6, ax7, ax8]:
    ax.set_xlabel('Latitude')

cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.1, pad=0.12, aspect=60, shrink=0.5)
cb1.set_label('Zonal wind speed, m/s')

plt.savefig(plot_dir + 'monsoon_psi_u.pdf', format='pdf')
plt.close()


