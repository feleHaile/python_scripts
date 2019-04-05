''' 
12/02/2019 Make comparable hms: half_land, half_10_land, half_10_land_tibet, and CMAP
'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh

# Set plotting directory
plot_dir = '/scratch/rg419/plots/asymmetry_paper/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


def pick_lons(data, lonin):
    #Find index range covering specified longitudes
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    return lons


def isca_monsoon_hm(run, ax, lonin, title):
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    #levels=np.arange(0.,15.,2.)
    levels=np.arange(2.,19.,2.)
    
    tickspace = [1, 19, 37, 55]
    labels = ['1st Jan', '1st Apr', '1st Jul', '1st Oct']
    
    lons = pick_lons(data,lonin)
    f1 = (data.precipitation*86400.).sel(lon=lons).mean('lon').plot.contourf(ax=ax, x='xofyear', y='lat', levels = levels, add_labels=False, extend='both', cmap='Blues', add_colorbar=False)
    ax.set_ylim(-45,45)
    ax.grid(True,linestyle=':')
    ax.set_yticks(np.arange(-30.,31.,30.))
    ax.set_title(title)
    ax.set_xticks(tickspace)
    ax.set_xticklabels(labels,rotation=25)
    
    data.close()
    return f1


def cmap_monsoon_hm(ax, lonin, title):
    
    data = xr.open_dataset('/disca/share/rg419/CMAP_precip.pentad.mean.nc')
    data.coords['pentad'] = (('time'), np.tile(np.arange(1,74),38))
    data = data.groupby('pentad').mean('time')
    levels=np.arange(2.,19.,2.)
     
    tickspace = [1, 19, 37, 55]
    labels = ['1st Jan', '1st Apr', '30th Jun', '28th Sept']
    
    lons = pick_lons(data,lonin)

    f1 = data.precip.sel(lon=lons).mean('lon').plot.contourf(ax=ax, x='pentad', y='lat', levels = levels, add_labels=False, extend='both', cmap='Blues', add_colorbar=False)
    ax.set_ylim(-45,45)
    ax.grid(True,linestyle=':')
    ax.set_yticks(np.arange(-30.,31.,30.))
    ax.set_title(title)
    ax.set_xticks(tickspace)
    ax.set_xticklabels(labels,rotation=25)
    
    data.close()
    return f1
    


# Set figure parameters
rcParams['figure.figsize'] = 15, 9
rcParams['font.size'] = 14

# Start figure with 9 subplots
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3, sharey='row')
axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]
    
top_row_run = 'half_land'
#isca_monsoon_hm(top_row_run, ax1, lonin=[100,140], title='100-140E')
#isca_monsoon_hm(top_row_run, ax2, lonin=[140,180], title='140-180E')
#f1 = isca_monsoon_hm(top_row_run, ax3, lonin=[180,200], title='180-200E')

isca_monsoon_hm(top_row_run, ax1, lonin=[120,150], title='120-150E')
isca_monsoon_hm(top_row_run, ax2, lonin=[150,180], title='150-180E')
f1 = isca_monsoon_hm(top_row_run, ax3, lonin=[180,210], title='180-210E')

mid_row_run = 'half_nh_land_tibet_ol8'
#isca_monsoon_hm(mid_row_run, ax4, lonin=[100,140], title='100-140E')
#isca_monsoon_hm(mid_row_run, ax5, lonin=[140,180], title='140-180E')
#f2 = isca_monsoon_hm(mid_row_run, ax6, lonin=[180,200], title='180-200E')

isca_monsoon_hm(mid_row_run, ax4, lonin=[120,150], title='120-150E')
isca_monsoon_hm(mid_row_run, ax5, lonin=[150,180], title='150-180E')
f2 = isca_monsoon_hm(mid_row_run, ax6, lonin=[180,210], title='180-210E')

#cmap_monsoon_hm(ax7, lonin=[60,100], title='South Asia (60-100E)')
#cmap_monsoon_hm(ax8, lonin=[100,140], title='East Asia (100-140E)')
#f3 = cmap_monsoon_hm(ax9, lonin=[140,160], title='Western North Pacific (140-160E)')

cmap_monsoon_hm(ax7, lonin=[70,100], title='South Asia (70-100E)')
cmap_monsoon_hm(ax8, lonin=[100,130], title='East Asia (100-130E)')
f3 = cmap_monsoon_hm(ax9, lonin=[130,160], title='Western North Pacific (130-160E)')


ax1.set_ylabel('Latitude')
ax4.set_ylabel('Latitude')
ax7.set_ylabel('Latitude')

plt.subplots_adjust(left=0.07, right=0.97, top=0.95, bottom=0.1, hspace=0.5, wspace=0.2)
cb1=fig.colorbar(f1, ax=[ax1, ax2, ax3], use_gridspec=True, orientation = 'vertical',fraction=0.05, pad=0.03, aspect=15, shrink=1)
cb1.set_label('Precipitation, mm/day')
cb1=fig.colorbar(f2, ax=[ax4, ax5, ax6], use_gridspec=True, orientation = 'vertical',fraction=0.05, pad=0.03, aspect=15, shrink=1)
cb1.set_label('Precipitation, mm/day')
cb1=fig.colorbar(f3, ax=[ax7, ax8, ax9], use_gridspec=True, orientation = 'vertical',fraction=0.05, pad=0.03, aspect=15, shrink=1)
cb1.set_label('Precipitation, mm/day')

# Save as a pdf
plt.savefig(plot_dir + 'monsoon_hms_isca_cmap_option1.pdf', format='pdf')
plt.close()

