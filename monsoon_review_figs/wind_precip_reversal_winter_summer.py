'''5/09/2018 Use JRA-55 and GPCP monthly data to plot 1979-2016 wind and precip difference between MJJAS and NDJFM'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams

plot_dir = '/scratch/rg419/plots/monsoon_review_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
    
# Load in data
data_u = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/ucomp_monthly/atmos_monthly_together.nc')
data_v = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/vcomp_monthly/atmos_monthly_together.nc')
data_p = xr.open_dataset('/disca/share/rg419/CMAP_precip.mon.mean.nc')
#name_temp = '/scratch/rg419/obs_and_reanalysis/GPCP_monthly/gpcp_cdr_v23rB1_y%04d_m%02d.nc'
#names = [name_temp % (y, m) for y in range(1979,2017) for m in range(1,13)]
#data_p = xr.open_mfdataset(names)
land_mask='/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc'
land = xr.open_dataset(land_mask)

#print(data_u.lon)
#print(data_p.longitude)
#print(land.longitude)

#print(data_u.sel(time=slice('1979','2016')))
data_u = data_u.sel(time=slice('1979','2016'))
data_v = data_v.sel(time=slice('1979','2016'))
data_p = data_p.sel(time=slice('1979','2016'))

# Make climatologies
u_clim = data_u.groupby('time.month').mean('time')
v_clim = data_v.groupby('time.month').mean('time')
p_clim = data_p.groupby('time.month').mean('time')
print('means taken')

u_MJJAS = u_clim.sel(month=[5,6,7,8,9]).mean('month') 
v_MJJAS = v_clim.sel(month=[5,6,7,8,9]).mean('month') 
p_MJJAS = p_clim.sel(month=[5,6,7,8,9]).mean('month') 

u_NDJFM = u_clim.sel(month=[11,12,1,2,3]).mean('month')
v_NDJFM = v_clim.sel(month=[11,12,1,2,3]).mean('month')
p_NDJFM = p_clim.sel(month=[11,12,1,2,3]).mean('month')

# Calculate MJJAS and NDJFM difference
u_diff = u_clim.sel(month=[5,6,7,8,9]).mean('month') - u_clim.sel(month=[11,12,1,2,3]).mean('month')
v_diff = v_clim.sel(month=[5,6,7,8,9]).mean('month') - v_clim.sel(month=[11,12,1,2,3]).mean('month')
p_diff = p_clim.sel(month=[5,6,7,8,9]).mean('month') - p_clim.sel(month=[11,12,1,2,3]).mean('month')

# Define monsoon region using precip
monsoon_nh = (p_diff.precip > 2.) & (p_clim.precip.sel(month=[5,6,7,8,9]).sum('month')/p_clim.precip.sum('month') > 0.55) 
monsoon_sh = (p_diff.precip < -2.) & (p_clim.precip.sel(month=[11,12,1,2,3]).sum('month')/p_clim.precip.sum('month') > 0.55) 
lats_nh = [p_diff.lat[i] for i in range(len(p_diff.lat)) if p_diff.lat[i] > 0]
lats_sh = [p_diff.lat[i] for i in range(len(p_diff.lat)) if p_diff.lat[i] < 0]

#lats_tropics = [p_diff.lat[i] for i in range(len(p_diff.lat)) if p_diff.lat[i] > -30. and p_diff.lat[i] < 30.]

#itcz_aug= np.zeros(len(p_clim.nlon),)
#itcz_feb= np.zeros(len(p_clim.nlon),)

#for i in range(len(p_clim.nlon)):
    #print(p_clim.latitude[p_clim.precip.sel(month=8)[:,i] == p_clim.precip.sel(month=8)[:,i].max('nlat')])
#    itcz_aug[i] = p_clim.latitude[p_clim.precip.sel(month=8)[:,i] == p_clim.precip.sel(month=8, nlat=lats_tropics)[:,i].max('nlat')].values
#    itcz_feb[i] = p_clim.latitude[p_clim.precip.sel(month=2)[:,i] == p_clim.precip.sel(month=2, nlat=lats_tropics)[:,i].max('nlat')].values
    
#print(itcz_aug)

# Make plots

# Set figure parameters
rcParams['figure.figsize'] = 9, 15
rcParams['font.size'] = 18

ref_arrow=5
arrowdir='uv'

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex='col')

# Contour plot of precipitation difference
f1 = p_diff.precip.plot.contourf(ax=ax1, x='lon', y='lat', cmap='RdBu', levels=np.arange(-10.,10.1,1.), add_colorbar=False, extend='both', add_labels=False)
monsoon_nh.sel(lat=lats_nh).plot.contour(ax=ax1, x='lon', y='lat', levels=np.arange(-1.,2.,1.), colors='m', add_colorbar=False, linewidths=2, add_labels=False)
monsoon_sh.sel(lat=lats_sh).plot.contour(ax=ax1, x='lon', y='lat', levels=np.arange(-1.,2.,1.), colors='m', add_colorbar=False, linewidths=2, add_labels=False)
land.lsm[0,:,:].plot.contour(ax=ax1, x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', linewidths=2.)
#plt.plot(p_diff.longitude, itcz_aug, 'r')
#plt.plot(p_diff.longitude, itcz_feb, 'b')
b = ax1.quiver(u_diff.lon[::5], u_diff.lat[::2], u_diff.var33.sel(lev=85000.)[::2,::5], v_diff.var34.sel(lev=85000.)[::2,::5], angles=arrowdir, scale=200.)
ax1.set_title('MJJAS - NDJFM')
ax1.quiverkey(b, 5.,65., ref_arrow, str(ref_arrow) + ' m/s', coordinates='data', fontproperties={'weight': 'bold', 'size': 10}, color='k', labelcolor='k', labelsep=0.03, zorder=10)

# Contour plot of MJJAS precipitation
f1 = p_MJJAS.precip.plot.contourf(ax=ax2, x='lon', y='lat', cmap='RdBu', levels=np.arange(-10.,10.1,1.), add_colorbar=False, extend='both', add_labels=False)
monsoon_nh.sel(lat=lats_nh).plot.contour(ax=ax2, x='lon', y='lat', levels=np.arange(-1.,2.,1.), colors='m', add_colorbar=False, linewidths=2, add_labels=False)
monsoon_sh.sel(lat=lats_sh).plot.contour(ax=ax2, x='lon', y='lat', levels=np.arange(-1.,2.,1.), colors='m', add_colorbar=False, linewidths=2, add_labels=False)
land.lsm[0,:,:].plot.contour(ax=ax2, x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', linewidths=2.)
#plt.plot(p_diff.longitude, itcz_aug, 'r')
#plt.plot(p_diff.longitude, itcz_feb, 'b')
b = ax2.quiver(u_MJJAS.lon[::5], u_MJJAS.lat[::2], u_MJJAS.var33.sel(lev=85000.)[::2,::5], v_MJJAS.var34.sel(lev=85000.)[::2,::5], angles=arrowdir, scale=200.)
ax2.set_title('MJJAS')
ax2.quiverkey(b, 5.,65., ref_arrow, str(ref_arrow) + ' m/s', coordinates='data', fontproperties={'weight': 'bold', 'size': 10}, color='k', labelcolor='k', labelsep=0.03, zorder=10)

# Contour plot of NDJFM precipitation
f1 = p_NDJFM.precip.plot.contourf(ax=ax3, x='lon', y='lat', cmap='RdBu', levels=np.arange(-10.,10.1,1.), add_colorbar=False, extend='both', add_labels=False)
monsoon_nh.sel(lat=lats_nh).plot.contour(ax=ax3, x='lon', y='lat', levels=np.arange(-1.,2.,1.), colors='m', add_colorbar=False, linewidths=2, add_labels=False)
monsoon_sh.sel(lat=lats_sh).plot.contour(ax=ax3, x='lon', y='lat', levels=np.arange(-1.,2.,1.), colors='m', add_colorbar=False, linewidths=2, add_labels=False)
land.lsm[0,:,:].plot.contour(ax=ax3, x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', linewidths=2.)
#plt.plot(p_diff.longitude, itcz_aug, 'r')
#plt.plot(p_diff.longitude, itcz_feb, 'b')
b = ax3.quiver(u_NDJFM.lon[::5], u_NDJFM.lat[::2], u_NDJFM.var33.sel(lev=85000.)[::2,::5], v_NDJFM.var34.sel(lev=85000.)[::2,::5], angles=arrowdir, scale=200.)
ax3.set_title('NDJFM')
ax3.quiverkey(b, 5.,65., ref_arrow, str(ref_arrow) + ' m/s', coordinates='data', fontproperties={'weight': 'bold', 'size': 10}, color='k', labelcolor='k', labelsep=0.03, zorder=10)
#ax3.quiverkey(b, 45.,-90., ref_arrow, str(ref_arrow) + ' m/s', coordinates='data', fontproperties={'weight': 'bold', 'size': 10}, color='k', labelcolor='k', labelsep=0.03, zorder=10)

for ax in [ax1,ax2,ax3]:
    ax.grid(True,linestyle=':')
    ax.set_ylabel('Latitude')
    ax.set_ylim(-60.,60.)
    ax.set_xticks(np.arange(0.,360.,90.))
    ax.set_yticks(np.arange(-60.,61.,30.))

ax1.text(-50, 60, 'a)')
ax2.text(-50, 60, 'b)')
ax3.text(-50, 60, 'c)')

ax3.set_xlabel('Longitude')

plt.subplots_adjust(left=0.12, right=0.97, top=0.97, bottom=0.02, hspace=0.2, wspace=0.1)

cb1=fig.colorbar(f1, ax=[ax1,ax2,ax3], use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.07, aspect=30, shrink=0.5)
cb1.set_label('Precipitation, mm/day')

plt.savefig(plot_dir + 'wind_precip_reversal_winter_summer.pdf', format='pdf')
plt.close()

data_u.close()
data_v.close()
data_p.close()