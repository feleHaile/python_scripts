# Function to generate land for experiments used in asymmetry paper

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import os


def other_land_shape(nlat, nlon, lat_array, lon_array):
    
    idx = np.zeros((nlat,nlon), dtype=bool)
    idx = (0. <= lon_array) & ((lon_array-160)*1.25 < lat_array) & (0. <= lat_array) & (lat_array <= 90.)   
    return idx

def write_land(exp, filename='land', land_mode='square',boundaries=[[20.,60.,20.,60.]],continents=['all'],topo_mode='none',mountains=['all'],topo_gauss=[[40.,40.,20.,10.,3500.]],waterworld=False):

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
    
    #create dictionary for continents
    cont_dic = {'NA':0, 'SA':1, 'EA':2, 'AF':3, 'OZ':4, 'IN':5, 'SEA':6}

# Firstly determine the land set-up to be used    
    # 1) Set-up in which a square of land is included
    if land_mode=='square':
        for blist in boundaries:
            idx = (blist[0] <= lat_array) & (lat_array <= blist[1]) & (blist[2] <= lon_array) & (blist[3] > lon_array)
            land_array[idx] = 1.0        
                
    elif land_mode=='none':  
        land_array = np.zeros((nlat,nlon))
    
    else:
        idx = other_land_shape(nlat, nlon, lat_array, lon_array)
        land_array[idx] = 1.0
        
    
# Next produce a topography array
    if topo_mode == 'none':
        topo_array = np.zeros((nlat,nlon))
        
    elif topo_mode == 'sauliere2012':
        # Rockys from Sauliere 2012
        h_0 = 2670.
        central_lon = 247.5
        central_lat = 40.
        L_1 = 7.5
        L_2 = 20.
        gamma_1 = 42.
        gamma_2 = 42.
        delta_1 = ((lon_array - central_lon)*np.cos(np.radians(gamma_1)) + (lat_array - central_lat)*np.sin(np.radians(gamma_1)))/L_1
        delta_2 = (-(lon_array - central_lon)*np.sin(np.radians(gamma_2)) + (lat_array - central_lat)*np.cos(np.radians(gamma_2)))/L_2
        h_arr_rockys = h_0 * np.exp(-(delta_1**2. + delta_2**2.))
        idx_rockys = (h_arr_rockys / h_0 > 0.05) #s make sure exponentials are cut at some point - use the value from p70 of Brayshaw's thesis. 
        
        # Andes based on Brayshaw/Sauliere design(very approximate)
        h_0 = 4080.
        central_lon = 290.
        central_lat = -19.
        L_1 = 3.75
        L_2 = 20.
        gamma_1 = 8.
        gamma_2 = 8.
        delta_1 = ((lon_array - central_lon)*np.cos(np.radians(gamma_1)) + (lat_array - central_lat)*np.sin(np.radians(gamma_1)))/L_1
        delta_2 = (-(lon_array - central_lon)*np.sin(np.radians(gamma_2)) + (lat_array - central_lat)*np.cos(np.radians(gamma_2)))/L_2
        h_arr_andes = h_0 * np.exp(-(delta_1**2. + delta_2**2.))
        idx_andes = (h_arr_andes / h_0 > 0.05) #s make sure exponentials are cut at some point - use the value from p70 of Brayshaw's thesis. 
        
        # West African highlands based on Brayshaw/Sauliere design (very approximate)
        h_0 = 1530.
        central_lon = 36.
        central_lat = -6.
        L_1 = 7.5
        L_2 = 36.
        gamma_1 = -36.
        gamma_2 = -36.
        delta_1 = ((lon_array - central_lon)*np.cos(np.radians(gamma_1)) + (lat_array - central_lat)*np.sin(np.radians(gamma_1)))/L_1
        delta_2 = (-(lon_array - central_lon)*np.sin(np.radians(gamma_2)) + (lat_array - central_lat)*np.cos(np.radians(gamma_2)))/L_2
        h_arr_afhi = h_0 * np.exp(-(delta_1**2. + delta_2**2.))
        idx_afhi = (h_arr_afhi / h_0 > 0.05) #s make sure exponentials are cut at some point - use the value from p70 of Brayshaw's thesis.
        
        #Tibet from Sauliere 2012
        h_0 = 5700.
        central_lon = 130. #82.5 #130.
        central_lat = 28
        L_1 = 12.5
        L_2 = 12.5
        gamma_1 = -49.5
        gamma_2 = -18.
        delta_1 = ((lon_array - central_lon)*np.cos(np.radians(gamma_1)) + (lat_array - central_lat)*np.sin(np.radians(gamma_1)))/L_1
        delta_2 = (-(lon_array - central_lon)*np.sin(np.radians(gamma_2)) + (lat_array - central_lat)*np.cos(np.radians(gamma_2)))/L_2
        h_arr_tibet_no_amp = np.exp(-(delta_1**2.))*(1./delta_2)*np.exp(-0.5*(np.log(delta_2))**2.)
        maxval = np.nanmax(h_arr_tibet_no_amp) #For some reason my maximum value of h_arr_tibet_no_amp > 1. Renormalise so h_0 sets amplitude. 
        h_arr_tibet = (h_arr_tibet_no_amp/maxval)*h_0
        idx_tibet = (h_arr_tibet / h_0 > 0.05)

        if 'all' in mountains:
            topo_array[idx_rockys] = h_arr_rockys[idx_rockys]
            topo_array[idx_afhi] =  h_arr_afhi[idx_afhi]
            topo_array[idx_tibet] =  h_arr_tibet[idx_tibet]
            topo_array[idx_andes] =  h_arr_andes[idx_andes]
        elif 'rockys' in mountains:
            topo_array[idx_rockys] = h_arr_rockys[idx_rockys]
        elif 'andes' in mountains:
            topo_array[idx_andes] =  h_arr_tibet[idx_andes]
        elif 'afhi' in mountains:
            topo_array[idx_afhi] =  h_arr_tibet[idx_afhi]
        elif 'tibet' in mountains:
            topo_array[idx_tibet] =  h_arr_tibet[idx_tibet]
        else:
            print('No valid mountain options detected for Sauliere 2012 topography')

            
    elif topo_mode == 'gaussian':
        #Options to define simple Gaussian Mountain
        for hill in topo_gauss:
            central_lat = hill[0]
            central_lon = hill[1]
            L_1 = hill[2]
            L_2 = hill[3]
            gamma = hill[4]
            h_0 = hill[5]
            delta_1 = ((lon_array - central_lon)*np.cos(np.radians(gamma)) + (lat_array - central_lat)*np.sin(np.radians(gamma)))/L_1
            delta_2 = (-(lon_array - central_lon)*np.sin(np.radians(gamma)) + (lat_array - central_lat)*np.cos(np.radians(gamma)))/L_2
            h_arr = h_0 * np.exp(-(delta_1**2. + delta_2**2.))
            idx = (h_arr / h_0 > 0.05)
            topo_array[idx] = h_arr[idx]
            
            #central_lat = hill[0]
            #central_lon = hill[1]
            #radius_degrees = hill[2]
            #std_dev = hill[3]
            #height = hill[4]
            #rsqd_array = np.sqrt((lon_array - central_lon)**2.+(lat_array - central_lat)**2.)
            #generalise to ellipse - needs checking but may be useful later (RG)
            #ax_rot = 1. #gradient of new x axis
            #ax_rat = 2. #axis ratio a**2/b**2
            #rsqd_array = np.sqrt((lon_array - central_lon + ax_rot*(lat_array - central_lat))**2.+ ax_rat*(lat_array - central_lat - ax_rot*(lon_array - central_lon))**2.)*np.cos(np.arctan(ax_rot))
            #divide by factor of cos(atan(m)) to account for change in coords
            #idx = (rsqd_array < radius_degrees) 
            #topo_array[idx] = height* np.exp(-(rsqd_array[idx]**2.)/(2.*std_dev**2.))
        
    else:
        print('Invalid topography option given')
        

    if waterworld != True:      #Leave flexibility to allow aquamountains!
        idx = (land_array == 0.) & (topo_array != 0.)
        topo_array[idx] = 0. 


    #Write land and topography arrays to file
    topo_filename = '/scratch/rg419/Experiments/' + exp + '/input/' + filename + '.nc'
    print(topo_filename)
    topo_file = Dataset(topo_filename, 'w', format='NETCDF3_CLASSIC')
    lat = topo_file.createDimension('lat', nlat)
    lon = topo_file.createDimension('lon', nlon)
    latitudes = topo_file.createVariable('lat','f4',('lat',))
    longitudes = topo_file.createVariable('lon','f4',('lon',))
    topo_array_netcdf = topo_file.createVariable('zsurf','f4',('lat','lon',))
    land_array_netcdf = topo_file.createVariable('land_mask','f4',('lat','lon',))
    latitudes[:] = lats
    longitudes[:] = lons
    topo_array_netcdf[:] = topo_array
    land_array_netcdf[:] = land_array
    topo_file.close()
    print(('Output written to: ' + topo_filename))


    #Show configuration on screen to allow checking
    lon_0 = lons.mean()
    lat_0 = lats.mean()
    m = Basemap(lat_0=lat_0,lon_0=lon_0)
    xi, yi = m(lon_array, lat_array)
    plt.figure()
    if land_mode != 'none':
        m.contour(xi,yi,land_array)
    if topo_mode != 'none':
        cs = m.contourf(xi,yi,topo_array, cmap=plt.get_cmap('RdBu_r'))
        cb = plt.colorbar(cs, shrink=0.5, extend='both')
    plt.xticks(np.linspace(0,360,13))
    plt.yticks(np.linspace(-90,90,7))
    plt.show()



if __name__ == "__main__":
    
    write_land('asym_aquaplanets', filename='half_nh_land_tibet_diag', land_mode='other', topo_mode='sauliere2012', mountains=['tibet'])
    