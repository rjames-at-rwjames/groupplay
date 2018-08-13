# To plot all CMIP5 models
# Scatterplots - so far options include :
#


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
import dsets_mplot_28 as dset_mp
import MetBot.mast_dset_dict as mast_dict
import MetBot.dimensions_dict as dim_exdict
import MetBot.mytools as my
import MetBot.MetBlobs as blb
import MetBot.EventStats as stats
import MetBot.mynetcdf as mync
import MetBot.SynopticAnatomy as sy
import scipy

### Which plot?
# index
index1='AngolaLow' # AngolaLow or Froude
seas1='DJF' # see text files for options

index2='Froude'
seas2='JF'

### Get directories
bkdir=cwd+"/../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
thisdir=bkdir+"groupplay/"
group=True

figdir=thisdir+"scatter_AL_Froude/"
my.mkdir_p(figdir)

### Dsets
dsets='all'
dsetnames=['noaa','cmip5']
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

print "Looping datasets"
for d in range(ndset):
    dset=dsetnames[d]
    dcnt=str(d+1)
    print 'Running on '+dset
    print 'This is dset '+dcnt+' of '+ndstr+' in list'

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


        ### Get index 1
        ind1=0
        indtxt1 = index1+'_'+seas1+'.txt'
        print 'Getting '+index1+' indices for this model'
        print indtxt1
        with open(indtxt1) as f:
            for line in f:
                if name in line:
                    ind1 = line.split()[1]
                    print 'index1=' + str(ind1)

        if ind1!=0:

            xvals[cnt]=ind1

            ### Get index 2
            ind2 = 0
            indtxt2 = index2 + '_' + seas2 + '.txt'
            print 'Getting ' + index2 + ' indices for this model'
            print indtxt2
            with open(indtxt2) as f:
                for line in f:
                    if name in line:
                        ind2 = line.split()[1]
                        print 'index2=' + str(ind2)

            if ind2 != 0:
                yvals[cnt] = ind2

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


xlab=index1+' '+seas1
ylab=index2+' '+seas2

plt.xlabel(xlab)
if index1=='AngolaLow':
    plt.xticks([1460,1470,1480,1490,1500])
plt.ylabel(ylab)

box=ax.get_position()
ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
ax.legend(loc='center left',bbox_to_anchor=[1,0.5],fontsize='xx-small',markerscale=0.8,numpoints=1)

plt.title('r2 '+str(round(rsquared,2)))



figsuf=""

if group:
    figsuf=figsuf+'_grouped'

scatterfig=figdir+'/scatter.'+index1+'.seas_'+seas1+'.'+index2+'.seas_'+seas2+'.'\
           +figsuf+'.png'
print 'saving figure as '+scatterfig
plt.savefig(scatterfig,dpi=150)
plt.close()