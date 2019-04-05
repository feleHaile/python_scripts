import matplotlib.pyplot as plt
import numpy as np
import shtns
import os
from matplotlib.ticker import FormatStrFormatter
from netCDF4 import Dataset

blank = FormatStrFormatter('')

prefix = 'walker_schneider_test'
run = '-rot10th'
tnum = 42
dirpath = '%s/T%d' % (prefix+run,tnum)
fieldname = 'vor'
filename_out = fieldname+'_spectral'+run+'.eps'
if os.path.exists(filename_out):
    print '! WARNING: %s will be overwritten' % (filename_out)
title = os.path.basename(__file__)+': '+dirpath+': '+fieldname+' by $p$ / hPa (this plot: '+filename_out+')'

w = 350

width = 11.69
height = 8.27
plt.figure(figsize=(width,height))
rows = 6
cols = 5

nfile = Dataset(filename=dirpath+'/atmos_daily.nc')

#time = nfile.variables['time']
#phalf = nfile.variables['phalf']
pfull = nfile.variables['pfull']
#lat = nfile.variables['lat']
#latb = nfile.variables['latb']
#lon = nfile.variables['lon']
#slp = nfile.variables['slp']                 # time,        lat, lon
#ps = nfile.variables['ps']                   # time,        lat, lon
#zsurf = nfile.variables['zsurf']             #              lat, lon
#pres_half   = nfile.variables['pres_half']   # time, phalf, lat, lon
#height_half = nfile.variables['height_half'] # time, phalf, lat, lon
#pres_full   = nfile.variables['pres_full']   # time, pfull, lat, lon
#height_full = nfile.variables['height']      # time, pfull, lat, lon
field = nfile.variables[fieldname]      # time, pfull, lat, lon

#field_window = field[-w:,:,:,:]
field_state = field[-1,:,:,:]	# snapshot for now

nlon = field.shape[3]
nlat = field.shape[2]

mmax = min(tnum,(nlon - 1)//2)
print tnum,mmax
sh = shtns.sht(tnum,mmax)	# reduce 2nd argument here as needed
#nlm = sh.nlm	# size of combined lm index
sh.set_grid(nlat,nlon)              # build default grid (gauss grid, phi-contiguous)

#print np.sctypeDict.keys()
field_k_l_m = np.empty((field.shape[1],tnum+1,mmax+1),dtype='complex64')
for k in range(0,field.shape[1]):
    field_lm = sh.analys(field_state[k,:,:].astype('float64'))	# sh.spec_array implicit?
    for l in range(0,tnum+1):
        for m in range(0,mmax+1):
            field_k_l_m[k,l,m] = field_lm[sh.idx(l,m)] if m <= l else np.nan

uspec = np.abs(field_k_l_m)**2	# pfull, l, m
uspec[:,:,1:] *= 2	# weight factor for missing m<0 (reality condition); 'optional' if not thinking of as summing displayed components
ushape = uspec.shape
print 'field shape:', ushape
#L = ushape[1]//2
L = min(20,ushape[1])
M = min(10,ushape[2])
#M = ushape[2]//3 + 2	# 2 for signed x 3/2 to eliminate anti-aliasing
#maxis = range(-M,M)
#print maxis
#maxis  = range(0,M)
#maxisp = range(0,M+1)	# +1 for pcolor
logmin = np.log10(np.nanmin(uspec[uspec > 0.]))
logmax = np.log10(np.nanmax(uspec))
print logmin,logmax
logsum = np.log10(np.nansum(uspec))
logminint = np.floor(logmin).astype(int)
logmaxint = np.ceil (logmax).astype(int)
print logminint,logmaxint
decs = range(logminint,logmaxint+1)	# +1 to accommodate later diff
dist = [np.count_nonzero(np.log10(uspec) < n) for n in decs]
print 'percentage per log band:', zip(decs,np.diff(dist)*100./uspec.size)
if dist[-1] != uspec.size:
    print '! incomplete'
logmin = max(logmin,logmax-9)	# override
logmin = max(logmin,logmax-2)	# override: less than 1% of max power = black
logmin = max(logmin,logsum-3)	# override: less than 0.1% of total power = black
print 'colorbar limits:', logmin, logmax

uspec = uspec[:,0:L,0:M]
#loguspec = np.ma.masked_where(loguspec < logmin,loguspec)
if M > L:
    loguspec = np.log10(uspec)
    xlabel = '$m$'
    ylabel = '$\\ell$'
else:
    loguspec = np.log10(np.swapaxes(uspec,1,2))
    xlabel = '$\\ell$'
    ylabel = '$m$'

protect_final = 2
psize = rows*cols - protect_final
usize = ushape[0] - protect_final	# NB ushape[1-2] out-of-date
if usize > psize:
    #print '! too many levels; will decimate'
    if usize >= 2*psize:
        mod = np.ceil(np.float(usize)/psize).astype(int)
        cap = 1
    else:
        mod = np.floor(1./((np.float(usize)/psize) - 1)).astype(int) + 1
        cap = mod - 1
    print '! levels plotted: first %d of every %d, plus final %d' % (cap,mod,protect_final)
else:
    mod = usize
    cap = mod
sub = 0
for k in range(0,ushape[0]):
    if k < usize and k % mod >= cap:
        print '!', k, pfull[k]
## RECALCULATE decimation here?
        continue
    sub += 1
#if sub == 1:
    ax = plt.subplot(rows,cols,sub)
    ax.set_aspect('equal')
#    ax1 = ax
#else:
#    ax = plt.subplot(rows,cols,sub,sharex=ax1,sharey=ax1)
    ax.set_title('%.2f' % (pfull[k]))
    #ax.set_yticks(np.linspace(-90,90,7))
    if (sub - 1) // cols + 1 == rows:
        ax.set_xlabel(xlabel)
#else:
#    ax.xaxis.set_major_formatter(blank)
#    plt.xticks([],[])
#    plt.yticks([],[])
    if sub % cols == 1:
        ax.set_ylabel(ylabel)
    #else:
    #    ax.yaxis.set_major_formatter(blank)

    c = ax.pcolor(loguspec[k,:,:],vmin=logmin,vmax=logmax,cmap='inferno',edgecolors='none')
    #c = ax.contourf(maxis,lat[:],loguspec[k,:,:],levels=np.linspace(logmin,logmax,3*9+1),cmap='inferno',extend='both')
    #plt.colorbar(c)
    print k, pfull[k] #, loguspec[k,:,:].min(), loguspec[k,:,:].max()

plt.suptitle(title)

plt.tight_layout()
plt.subplots_adjust(top=0.9)

plt.savefig(filename_out)
