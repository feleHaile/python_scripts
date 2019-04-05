'''30/11/2018 Load COBE SST data and evaluate niño3.4 index'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
import datetime
import pandas as pd

# Set and make plotting directory
plot_dir = '/scratch/rg419/plots/onset_variability_new/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

# Load in individual datafiles and save to single file
#path_template = '/scratch/rg419/obs_and_reanalysis/cobe_sst/sst%04d%02d.grb'
#paths = [path_template % (y, m) for y in range(1943,2017) for m in range(1,13)] #2017
#data = xr.open_mfdataset(paths, engine='cfgrib', concat_dim='time')
#data.to_netcdf('/scratch/rg419/obs_and_reanalysis/cobe_sst/cobe_sst_1943_2016.nc')

# Load in data
data = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/cobe_sst/cobe_sst_1943_2016.nc')

# Weight SST by cos(lat) for area mean
coslat = np.cos(data.latitude * np.pi/180.)
sst_weighted = data.p3080 * coslat # Weight sst by coslat, to do area weighting

# Make a list of lats and lons in the Nino 3.4 region
lats = [data.latitude[i].values for i in range(len(data.latitude)) if data.latitude[i] >= -5. and data.latitude[i] <= 5.]
lons = [data.longitude[i].values for i in range(len(data.longitude)) if data.longitude[i] >= 190. and data.longitude[i] <= 240.]

# Calculate area mean SST
nino_34_sst = sst_weighted.sel(longitude=lons).sel(latitude=lats).sum(('latitude','longitude')) / (coslat.sel(latitude=lats).sum('latitude') * len(lons))

# Take 30 year rolling mean centred on each year
nino_34_rm = nino_34_sst.rolling(time=360, center=True).mean()
# Cut blank years at the start where rolling mean was not defined
nino_34_rm = nino_34_rm[180:]
# For later part of data where rolling mean is not defined, use the mean over the last 30 years of the data set
nino_34_rm[-179:] = nino_34_sst[-360:].mean('time')

# Evaluate nino34 index - subtract rollng mean from monthly data
nino_34 = nino_34_sst[180:] - nino_34_rm
nino_34 = nino_34.drop(['depthBelowSea','step','valid_time'])
nino_34.name = 'nino_34'
# Take a year long rolling mean to smooth to make nice plots, if wanted. NB for monthly index need unsmoothed data
#nino_34 = nino_34.rolling(time=12, center=True).mean()

#nino_34 = nino_34_sst[180:] - nino_34_sst[180:].mean('time')
nino_34.to_netcdf('/scratch/rg419/obs_and_reanalysis/nino_34.nc')

# Get march nino 3.4 index: Can be done by either indexing using a datetime object, or adding a new coord. The former is simpler but the latter more general
nino_34_jan = nino_34.sel(time=pd.to_datetime(['%04d-01-15' % y for y in range(1958,2017)]))
nino_34_march = nino_34.sel(time=pd.to_datetime(['%04d-03-15' % y for y in range(1958,2017)]))
#nino_34_march = nino_34.assign_coords(month = ('time', nino_34['time.month'])).set_index(time=['time','month']).sel(month=3)

def nino_plot(nino, title=''):
    nino.plot(x='time', color='k')
    plt.plot([nino.time.values[0],nino.time.values[-1]],[0,0], color='k')
    plt.fill_between(nino.time.values, 0., nino, where=nino>=0., color='C1')
    plt.fill_between(nino.time.values, 0., nino, where=nino<=-0., color='C0')
    
    plt.plot([nino.time.values[0],nino.time.values[-1]],[0.5,0.5], color='r', linestyle='--')
    plt.plot([nino.time.values[0],nino.time.values[-1]],[-0.5,-0.5], color='b', linestyle='--')
    plt.ylabel('Anomaly in Degrees C')
    plt.xlabel('Year')
    plt.title('SST Anomaly in Nino 3.4 Region')
    plt.savefig(plot_dir + 'nino_34' + title + '.pdf', format='pdf')
    plt.close()

nino_plot(nino_34)
nino_plot(nino_34_jan, '_jan')
nino_plot(nino_34_march, '_march')
