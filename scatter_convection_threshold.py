# To plot all CMIP5 models or UM models
# Scatterplots

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma

cwd=os.getcwd()
sys.path.append(cwd+'/..')
sys.path.append(cwd+'/../MetBot')
sys.path.append(cwd+'/../RTools')
sys.path.append(cwd+'/../quicks')
import MetBot.dset_dict as dsetdict
import MetBot.mast_dset_dict as mast_dict
import MetBot.dimensions_dict as dim_exdict
import MetBot.mytools as my
import MetBot.MetBlobs as blb
import MetBot.EventStats as stats
import MetBot.mynetcdf as mync
import MetBot.SynopticAnatomy as sy
import scipy

# Which dataset?

whichd='CMIP5' # UM or CMIP5

if whichd=='UM':
    import dsets_mplot_um as dset_mp
elif whichd=='CMIP5':
    import dsets_mplot_28 as dset_mp
group=True # works for CMIP5 only

### What are we plotting?

# season
seas='DJF'

# x axis - convection
globv1='omega' # olr or omega
levsel1=True
if levsel1:
    choosel1=['500'] # can add a list
else:
    choosel1=['1']
sub_x='subt'
l=0

# y axis - threshold
thname='actual'



### Get directories
bkdir=cwd+"/../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
thisdir=bkdir+"groupplay/"
refkey='0'

figdir=thisdir+"scatter_convection_thresh/"
my.mkdir_p(figdir)

### Dsets
dsets='all'
if whichd=='CMIP5':
    dsetnames=['noaa','cmip5']
elif whichd=='UM':
    dsetnames = ['noaa', 'cmip5', 'um']
ndset=len(dsetnames)
ndstr=str(ndset)

### Count total number of models
nm_dset=np.zeros(ndset)
for d in range(ndset):
    dset=dsetnames[d]
    nmod=len(dset_mp.dset_deets[dset])
    nm_dset[d]=nmod
nallmod=np.sum(nm_dset)
nallmod=int(nallmod)
print nallmod


print "Setting up plot..."
plt.figure(figsize=[6,5])
ax=plt.subplot(111)

### Get arrays for output
modnames=[]

### colours
if not group:
    cols=['b','g','r','c','m','gold','k',\
        'b','g','r','c','m','gold','k',\
        'b','g','r','c','m','gold','k',\
        'b','g','r','c','m','gold','k']
    markers=["o","o","o","o","o","o","o",\
        "^","^","^","^","^","^","^",\
        "*","*","*","*","*","*","*",\
        "d","d","d","d","d","d","d"]
elif group:
    grcls=['fuchsia','b','r','blueviolet','springgreen','gold','darkorange']
    grcnt=np.zeros(7,dtype=np.int8)
    grmrs=["o","^","*","d","+","v","h"]

## other plotting set up
siz=np.full((28),5)
siz[0]=8
xvals=np.ma.zeros(nallmod,dtype=np.float32)
yvals=np.ma.zeros(nallmod,dtype=np.float32)
cnt = 0

### Get season info
if seas == 'NDJFM':
    mons = [1, 2, 3, 11, 12]
    mon1 = 11
    mon2 = 3
    nmon = 5
elif seas == 'DJF':
    mons = [1, 2, 12]
    mon1 = 12
    mon2 = 2
    nmon = 3
elif seas == 'JF':
    mons = [1, 2]
    mon1 = 1
    mon2 = 2
    nmon = 2


print "Looping datasets"
for d in range(ndset):
    dset=dsetnames[d]
    dcnt=str(d+1)
    print 'Running on '+dset
    print 'This is dset '+dcnt+' of '+ndstr+' in list'

    if dset != 'cmip5': levc = int(choosel1[l])
    else: levc = int(choosel1[l]) * 100

    ### Models
    mods = 'all'
    nmod = len(dset_mp.dset_deets[dset])
    mnames = list(dset_mp.dset_deets[dset])
    nmstr = str(nmod)

    for mo in range(nmod):
        name = mnames[mo]
        mcnt = str(mo + 1)
        print 'Running on ' + name
        print 'This is model ' + mcnt + ' of ' + nmstr + ' in list'

        ### if group get group
        if group:
            groupdct = dset_mp.dset_deets[dset][name]
            thisgroup = int(groupdct['group'])
            grcl = grcls[thisgroup - 1]
            grmr=grmrs[grcnt[thisgroup-1]]
            grcnt[thisgroup-1]+=1


        ### Convection information for x axis
        # Switch variable if NOAA
        if dset == 'noaa' and globv1 != 'olr':
            if globv1 == 'pr':
                ds4noaa = 'trmm'
                mod4noaa = 'trmm_3b42v7'
            else:
                ds4noaa = 'ncep'
                mod4noaa = 'ncep2'
            dset2 = ds4noaa
            name2 = mod4noaa
        else:
            dset2 = dset
            name2 = name

        # Get info
        moddct = dsetdict.dset_deets[dset2][name2]
        moddct1 = dset_mp.dset_deets[dset2][name2]
        vnamedict = globv1 + 'name'
        mastdct = mast_dict.mast_dset_deets[dset2]
        varstr = mastdct[vnamedict]
        dimdict = dim_exdict.dim_deets[globv1][dset2]
        latname = dimdict[1]
        lonname = dimdict[2]
        if dset2 == 'um':
            ys = moddct1['climyr']
        else:
            if globv1 != 'omega' and globv1 != 'q' and globv1 != 'gpth':
                ys = moddct['yrfname']
            else:
                if name2 == "MIROC5":
                    if globv1 == 'q':
                        ys = moddct['fullrun']
                    elif globv1 == 'omega' or globv1 == 'gpth':
                        ys = '1950_2009'
                    else:
                        print 'variable ' + globv1 + ' has unclear yearname for ' + name2
                else:
                    ys = moddct['fullrun']

        # Open ltmonmean file
        meanfile = bkdir + 'metbot_multi_dset/' + dset2 + '/' + name2 + '/' \
                   + name2 + '.' + globv1 + '.mon.mean.' + ys + '.nc'

        print 'Attempting to open ' + meanfile

        if os.path.exists(meanfile):

            if levsel1:
                ncout = mync.open_multi(meanfile, globv1, name2, \
                                        dataset=dset2, subs=sub_x, levsel=levc)
            else:
                ncout = mync.open_multi(meanfile, globv1, name2, \
                                        dataset=dset2, subs=sub_x)
            print '...file opened'
            ndim = len(ncout)
            if ndim == 5:
                meandata, time, lat, lon, dtime = ncout
            elif ndim == 6:
                meandata, time, lat, lon, lev, dtime = ncout
                meandata = np.squeeze(meandata)
            else:
                print 'Check number of dims in ncfile'
            dtime[:, 3] = 0

            # Remove duplicate timesteps
            print 'Checking for duplicate timesteps'
            tmp = np.ascontiguousarray(dtime).view(
                np.dtype((np.void, dtime.dtype.itemsize * dtime.shape[1])))
            _, idx = np.unique(tmp, return_index=True)
            dtime = dtime[idx]
            meandata = meandata[idx, :, :]

            nlat = len(lat)
            nlon = len(lon)

            # Select seasons and get mean
            thesemons = np.zeros((nmon, nlat, nlon), dtype=np.float32)
            for zz in range(len(mons)):
                thesemons[zz, :, :] = meandata[mons[zz] - 1, :, :]
            seasmean = np.nanmean(thesemons, 0)

            # Get basic regional mean
            reg_mean=np.nanmean(seasmean)

            xvals[cnt]=reg_mean


            ### TTT thresh for y axis
            ### Get threshold for TTTs
            threshtxt = botdir + 'thresholds.fmin.all_dset.txt'
            print 'Getting threshold for this model'
            print threshtxt
            with open(threshtxt) as f:
                for line in f:
                    if dset + '\t' + name in line:
                        thresh = line.split()[2]
                        print 'thresh=' + str(thresh)

            thresh = int(thresh)
            thisthresh = thresh
            thre_str = str(int(thisthresh))


            yvals[cnt]=thisthresh

            print name
            print xvals[cnt]
            print yvals[cnt]
            if group:
                colour = grcl
                mk = grmr
            else:
                colour = cols[cnt]
                mk = markers[cnt]

            ax.plot(xvals[cnt], yvals[cnt], marker=mk, \
                    color=colour, label=name, markeredgecolor=colour, markersize=siz[cnt], linestyle='None')

        else:
            xvals[cnt] = ma.masked
            yvals[cnt] = ma.masked

        cnt += 1
        modnames.append(name)

# Set up plot
print xvals
print yvals
grad, inter, r_value, p_value, std_err = scipy.stats.mstats.linregress(xvals, yvals)
rsquared = r_value ** 2
if rsquared > 0.4:
    ax.plot(xvals, (grad * xvals + inter), '-', color='k')


xlab=globv1+'_'+sub_x+'_'+seas

ylab = 'OLR threshold'

plt.xlabel(xlab)
plt.ylabel(ylab)

box=ax.get_position()
ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
ax.legend(loc='center left',bbox_to_anchor=[1,0.5],fontsize='xx-small',markerscale=0.8,numpoints=1)

plt.title('r2 '+str(round(rsquared,2)))

figsuf=""

if group:
    figsuf=figsuf+'_grouped'

scatterfig=figdir+'/scatter.'+globv1+'.'+sub_x+'.seas_'+seas+'.OLRthreshold.'\
           +thname+'.png'
print 'saving figure as '+scatterfig
plt.savefig(scatterfig,dpi=150)
plt.close()