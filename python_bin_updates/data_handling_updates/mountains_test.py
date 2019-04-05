# Test gaussian mountains

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import os

def mountain_test(hill = [90., 0., 10., 50., 60., 40., 3500.], i=1):

# Common features of set-ups
    # specify resolution
    t_res = 42
    #read in grid from approriate file
    GFDL_BASE = os.environ['GFDL_BASE']
    resolution_file = Dataset(GFDL_BASE + 'src/extra/python/scripts/gfdl_grid_files/t'+str(t_res)+'.nc', 'r', format='NETCDF3_CLASSIC')
    lons = resolution_file.variables['lon'][:]
    lats = resolution_file.variables['lat'][:]
    lonb = resolution_file.variables['lonb'][:]
    latb = resolution_file.variables['latb'][:]
    nlon=lons.shape[0]
    nlat=lats.shape[0]
    topo_array = np.zeros((nlat,nlon))
    land_array = np.zeros((nlat,nlon))

    #make 2d arrays of latitude and longitude
    lon_array, lat_array = np.meshgrid(lons,lats)
    lonb_array, latb_array = np.meshgrid(lonb,latb)
    
    central_lon = hill[0]
    central_lat = hill[1]
    L_1 = hill[2]
    L_2 = hill[3]
    gamma_1 = hill[4]
    gamma_2 = hill[5]
    h_0 = hill[6]
    
    delta_1 = ((lon_array - central_lon)*np.cos(np.radians(gamma_1)) + (lat_array - central_lat)*np.sin(np.radians(gamma_1)))/L_1
    delta_2 = (-(lon_array - central_lon)*np.sin(np.radians(gamma_2)) + (lat_array - central_lat)*np.cos(np.radians(gamma_2)))/L_2
    h_arr = h_0 * np.exp(-(delta_1**2. + delta_2**2.))
    idx = (h_arr / h_0 > 0.05) #s make sure exponentials are cut at some point - use the value from p70 of Brayshaw's thesis. 
    
    topo_array[idx] = h_arr[idx]
    
    #elif topo_mode == 'gaussian':
        #Options to define simple Gaussian Mountain
    #    for hill in topo_gauss:
    #        central_lat = hill[0]
    #        central_lon = hill[1]
    #        radius_degrees = hill[2]
    #        std_dev = hill[3]
    #        height = hill[4]
    #        rsqd_array = np.sqrt((lon_array - central_lon)**2.+(lat_array - central_lat)**2.)
    #        idx = (rsqd_array < radius_degrees) 
    #        topo_array[idx] = height* np.exp(-(rsqd_array[idx]**2.)/(2.*std_dev**2.))

    
    #Show configuration on screen to allow checking
    plt.figure(i)
    lon_0 = lons.mean()
    lat_0 = lats.mean()
    m = Basemap(lat_0=lat_0,lon_0=lon_0)
    xi, yi = m(lon_array, lat_array)
    cs = m.contourf(xi,yi,topo_array, cmap=plt.get_cmap('RdBu_r'))
    cb = plt.colorbar(cs, shrink=0.5, extend='both')
    plt.xticks(np.linspace(0,360,13))
    plt.yticks(np.linspace(-90,90,7))


if __name__ == "__main__":
    
    mountain_test(hill = [90., 0., 20., 20., 120., 30., 3500.], i=1)
    mountain_test(hill = [90., 0., 20., 20., 30., 120., 3500.], i=2)
    
    plt.show()
    