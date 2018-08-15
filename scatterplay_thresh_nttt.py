# Scatterplot wrapper
#   to plot number of events for each model versus threshold
#   we had done this before but now I am repeating with groups


# .....dset: noaa, um, cmip5
# .....name: model names will be taken from dset dictionary
# .....directory: here ../../CTdata/metbot_multi_dset/$dset/

import numpy as np
import matplotlib.pyplot as plt
import sys,os
cwd=os.getcwd()
sys.path.append(cwd+'/..')
sys.path.append(cwd+'/../MetBot')
sys.path.append(cwd+'/../RTools')
sys.path.append(cwd+'/../quicks')
import MetBot.SynopticAnatomy as sy
import MetBot.EventStats as stats
import MetBot.SynopticPlot as syp
import MetBot.AdvancedPlots as ap
import MetBot.RainStats as rs
import MetBot.MetBlobs as blb
import MetBot.mytools as my
import MetBot.mynetcdf as mync
import glob, socket, os
import mpl_toolkits.basemap as bm
import MetBot.dset_dict as dsetdict
import MetBot.find_saddle as fs
import scipy
import dsets_mplot_28 as dset_mp


### Running options
testyear=False  # plot based on 1 year of test data
threshtest=False # Option to run on thresholds + and - 5Wm2 as a test
group=True
wh_count = 'blob'  # blob or event

### Directory
bkdir=cwd+"/../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
thisdir=bkdir+"groupplay/"

figdir=thisdir+"scatter_thresh_nttt/"
my.mkdir_p(figdir)


### Loop threshs
if threshtest:
    thnames=['lower','actual','upper']
else:
    thnames=['actual']

nthresh=len(thnames)

### Loop domains
doms = ['All', 'Cont', 'Mada']
ndoms = len(doms)

for do in range(ndoms):

    for t in range(nthresh):

        ### Multi dset?
        dsets='all'     # "all" or "spec" to choose specific dset(s)
        dsetnames = ['noaa', 'cmip5']
        ndset = len(dsetnames)
        dsetstr=('_'.join(dsetnames))+'_'+str(ndset)
        print 'Running on datasets:'
        print dsetnames

        ### Count total number of models
        nm_dset=np.zeros(ndset)
        for d in range(ndset):
            dset = dsetnames[d]
            nmod = len(dset_mp.dset_deets[dset])
            nm_dset[d]=nmod
        nallmod=np.sum(nm_dset)
        nallmod=int(nallmod)
        print 'Total number of models = '+str(nallmod)

        ### colours
        if not group:
            cols = ['b', 'g', 'r', 'c', 'm', 'gold', 'k', \
                    'b', 'g', 'r', 'c', 'm', 'gold', 'k', \
                    'b', 'g', 'r', 'c', 'm', 'gold', 'k', \
                    'b', 'g', 'r', 'c', 'm', 'gold', 'k']
            markers = ["o", "o", "o", "o", "o", "o", "o", \
                       "^", "^", "^", "^", "^", "^", "^", \
                       "*", "*", "*", "*", "*", "*", "*", \
                       "d", "d", "d", "d", "d", "d", "d"]
        elif group:
            grcls = ['fuchsia', 'b', 'r', 'blueviolet', 'springgreen', 'gold', 'darkorange']
            grcnt = np.zeros(7, dtype=np.int8)
            grmrs = ["o", "^", "*", "d", "+", "v", "h"]

        ## other plotting set up
        siz = np.full((28), 5)
        siz[0] = 8
        xvals = np.ma.zeros(nallmod, dtype=np.float32)
        yvals = np.ma.zeros(nallmod, dtype=np.float32)
        cnt = 0

        ### Open figure
        plt.figure(figsize=[7, 5])
        ax = plt.subplot(111)


        ### Open arrays for results
        modnm=["" for x in range(nallmod)] # creates a list of strings for modnames

        ### Loop dsets and models
        for d in range(ndset):
            dset=dsetnames[d]
            nmod=len(dset_mp.dset_deets[dset])
            mnames=list(dset_mp.dset_deets[dset])
            print 'Looping through models'


            print mnames

            for mo in range(nmod):
                name=mnames[mo]
                mcnt = str(mo + 1)

                ### if group get group
                if group:
                    groupdct = dset_mp.dset_deets[dset][name]
                    thisgroup = int(groupdct['group'])
                    grcl = grcls[thisgroup - 1]
                    grmr = grmrs[grcnt[thisgroup - 1]]
                    grcnt[thisgroup - 1] += 1


                ### Get threshold
                if testyear:
                    threshtxt = botdir + 'thresholds.fmin.'+dset+'.test.txt'
                else:
                    threshtxt=botdir+'thresholds.fmin.all_dset.txt'
                print threshtxt
                with open(threshtxt) as f:
                    for line in f:
                        if dset+'\t'+name in line:
                            thresh = line.split()[2]
                            print 'thresh='+str(thresh)

                thresh = int(thresh)

                if thnames[t]=='actual':
                    thisthresh=thresh
                if thnames[t]=='lower':
                    thisthresh=thresh - 5
                if thnames[t]=='upper':
                    thisthresh=thresh + 5

                thre_str = str(int(thisthresh))
                xvals[cnt] = thisthresh



                # Find TTT data
                outsuf=botdir+dset+"/"+name+"/"+name+"_"

                mbsfile = outsuf + thre_str + '_' + dset + "-olr-0-0.mbs"
                syfile = outsuf + thre_str + '_' + dset + '-OLR.synop'


                # using events
                if wh_count=='event':

                    ###  Open synop file
                    s = sy.SynopticEvents((),[syfile],COL=False)

                    ### Count number of events
                    ks = s.events.keys() # all events
                    kw, ke = stats.spatialsubset(s,False,cutlon=40.) # splitting tracks west and east of 40E

                    ### Put n events into array
                    if doms[do]=='All':
                        yvals[cnt]=len(ks)
                    elif doms[do]=='Cont':
                        yvals[cnt]=len(kw)
                    elif doms[do]=='Mada':
                        yvals[cnt]=len(ke)


                elif wh_count=='blob':

                    refmbs, refmbt, refch = blb.mbopen(mbsfile)

                    blob_edts = []
                    blob_edts_cont =[]
                    blob_edts_mada = []
                    for b in range(len(refmbt)):
                        date = refmbt[b]
                        mon = int(date[1])
                        cX = refmbs[b, 3]
                        cY = refmbs[b, 4]

                        blob_edts.append(date)

                        if cX < 40.0:

                            blob_edts_cont.append(date)

                        elif cX >= 40.0:

                            blob_edts_mada.append(date)

                    blob_edts=np.asarray(blob_edts)
                    blob_edts_cont=np.asarray(blob_edts_cont)
                    blob_edts_mada=np.asarray(blob_edts_mada)

                    ### Put n events into array
                    if doms[do]=='All':
                        yvals[cnt]=len(blob_edts)
                    elif doms[do]=='Cont':
                        yvals[cnt]=len(blob_edts_cont)
                    elif doms[do]=='Mada':
                        yvals[cnt]=len(blob_edts_mada)


                ### Plot
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


                ### Put name into string list
                modnm[cnt]=dset+"_"+name
                cnt+=1

                print 'cnt'
                print cnt

        print xvals
        print yvals
        grad, inter, r_value, p_value, std_err = scipy.stats.mstats.linregress(xvals, yvals)
        rsquared=r_value**2

        if rsquared > 0.4:
            ax.plot(xvals, (grad * xvals + inter), '-', color='k')

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=[1, 0.5], fontsize='xx-small', markerscale=0.8, numpoints=1)

        plt.title('r2 ' + str(round(rsquared, 2)))

        plt.xlabel('OLR threshold')
        plt.ylabel('Number of TTT events')

        scatterfig=figdir+'/scatter_threshold_nTTT.thresh_'+thnames[t]+'.'+dsetstr+'.'+doms[do]+'.png'
        plt.savefig(scatterfig,dpi=150)
        plt.close()
