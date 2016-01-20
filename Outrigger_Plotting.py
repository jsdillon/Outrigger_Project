import numpy as np
import healpy as hp
import math
import sys,os,os.path
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from MAPS21cm.Specifications import Specifications
from MAPS21cm import Geometry
import cPickle as pickle
import scipy.constants as const
#from mpldatacursor import datacursor
from os import listdir
from os.path import isfile, join
import pylab
plt.close('all')

scriptDirectory = os.path.dirname(os.path.abspath(__file__))
allTests = pickle.load(open(scriptDirectory + "/Results/allTests.p","rb"))
#%%

def plotArray(antennaPositions):
    fig = plt.gcf()
    xLim = 1.2 * np.max(antennaPositions)
    yLim = xLim
    plt.xlim(-xLim,xLim)
    plt.ylim(-yLim,yLim)
    for antennaPos in antennaPositions:
        fig.gca().add_artist(plt.Circle((antennaPos[0],antennaPos[1]),7,fc='0.3',ec='k'))

 
def SortedEigensystem(matrix):
    """Returns the eigensystem of the input matrix where eigenvalues and eigenvectors are sorted by descending absolute value."""
    evals,evecs = np.linalg.eig(matrix)
    indices = np.argsort(np.abs(evals))[::-1]   
    return evals[indices], evecs[:,indices]


def dishScatter(fig, ax, xPos, yPos, cVals, colormap, radii=7.0, cLim = None, blackOutlines = 0): 
    if blackOutlines: radii -= .5
    if not hasattr(radii,'len'): radii = np.ones(len(xPos))*radii    
    if cLim is None:     
        underlyingScatter = plt.scatter(xPos, yPos, c=cVals, s=0, cmap=colormap.name)
        cVals = (cVals - np.min(cVals))/(np.max(cVals) - np.min(cVals))
    else:
        underlyingScatter = plt.scatter(xPos, yPos, c=cVals, s=0, cmap=colormap.name, vmin = cLim[0], vmax = cLim[1])
        cVals = 1.0*(cVals - cLim[0])/(cLim[1]-cLim[0])
    for x,y,r,c in zip(xPos, yPos, radii, colormap(cVals)): 
        ax.add_patch(plt.Circle((x,y), r, fc=c, ec='k',lw=blackOutlines))
    return underlyingScatter


#%% Load results of calcuations
allResults = {}
for test in allTests:
    #if test['title'].find('Long Integration') < 0 or test['title'].find('High-Res') < 0:
    if test['title'].find('Nine-Way') < 0:
        allResults[test['title']] = {}
        resultsDir =  scriptDirectory + '/Results/' + test['folder']
        onlyFiles = [ f for f in listdir(resultsDir) if isfile(join(resultsDir,f)) ]
        print 'Now loading ' + test['title']
        for f in onlyFiles:      
            if True:#f != 'mapNoiseCovariance.npy' and f != 'AtransNinvA.npy' and f != 'noiseEvecs.npy':# and f != 'PSF.npy':
                if f.split('.')[1] == 'npy': allResults[test['title']][f.split('.')[0]] = np.load(resultsDir + f)
                if f.split('.')[1] == 'p': allResults[test['title']][f.split('.')[0]] = pickle.load(open(resultsDir + f,'rb'))
                if f.split('.')[1] == 'txt': allResults[test['title']][f.split('.')[0]] = np.loadtxt(resultsDir + f)
#        print allResults[test['title']]['specifications'].nAntennas 
        

#%% Inrigger Calibratability

def CalibratabilityComparison(tests,figNum):
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12,10), num=figNum)
    for ax,test in zip(axes.flatten(), tests):
        s = allResults[test]['specifications']
        hexScatter = dishScatter(fig, ax, s.antennaPositions[:,0], s.antennaPositions[:,1], np.sqrt(allResults[test]['gainVariances']), pylab.cm.jet, 7.0, cLim=[1.0,3.0])
        ax.set_title(str(int(allResults[test]['nZeroEVs'][0])) + ' missing gain modes.\n' +  str(int(allResults[test]['nZeroEVs'][1])) + ' missing phase modes.')
        ax.axis('equal')    
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    cax = plt.axes([0.85, 0.1, 0.05, 0.8])
    fig.colorbar(hexScatter, cax=cax, cmap = pylab.cm.jet.name)


inriggerTests = ['HERA-331 + Fiducial Inriggers Full Sky',
         'HERA-331 + 6 Corner Inrigger Pairs Full Sky',
         'HERA-331 Split Core Full Sky']#,
#         'HERA-331 Nine-Way Split Core Full Sky']
CalibratabilityComparison(inriggerTests,2)


fullArrayTests = ['HERA-331 + Fidicuial Outriggers + Fidicuial Inriggers High-Res',
         'HERA-331 + Redundant Outriggers + 6 Corner Inrigger Pairs High-Res',
         'HERA-331 Split Core + Split Outriggers High-Res']#,
#         'HERA-331 Nine-Way Split Core + Split Outriggers High-Res']
CalibratabilityComparison(fullArrayTests,3)      





#%% High res mapmaking

highResTests = ['HERA-331 High-Res',
                'HERA-331 + Fidicuial Outriggers + Fidicuial Inriggers High-Res',
                'HERA-331 + Redundant Outriggers + 6 Corner Inrigger Pairs High-Res',
                'HERA-331 Split Core + Split Outriggers High-Res']#,
#                'HERA-331 Nine-Way Split Core + Split Outriggers High-Res']
plt.figure(5, figsize=(12,8)); plt.clf()
for test in highResTests:
    s = allResults[test]['specifications']        
    plt.semilogy(np.real(allResults[test]['noiseEvals'])/ s.nAntennas**2 * 331**2)
plt.xlabel('Eigenvalue Number',fontsize=16)
plt.ylabel('Eigenvalue of $\mathbf{A}^\dagger\mathbf{N}^{-1}\mathbf{A}$',fontsize=16)
plt.legend(highResTests, loc='upper right', fontsize=11)
plt.ylim([1e-6, 1e4])
plt.title('High Res Mapmaking')


#%% Full sky PSFs (This turns out to be rather unenlightening)

if False:
    fullSkyTests = ['HERA-331 Full Sky',
                    'HERA-331 Full Sky Long Integration',                
                    'HERA-331 + Fiducial Inriggers Full Sky',
                    'HERA-331 + Fiducial Inriggers Full Sky Long Integration',                
                    'HERA-331 + 6 Corner Inrigger Pairs Full Sky',
                    'HERA-331 + 6 Corner Inrigger Pairs Full Sky Long Integration',
                    'HERA-331 Split Core Full Sky',
                    'HERA-331 Split Core Full Sky Long Integration']
    
    
    fullSkyCoords = allResults['HERA-331 Full Sky']['coords']
    HERAFullSkyPSF = allResults['HERA-331 Full Sky']['PSF'][fullSkyCoords.newPSFIndices,fullSkyCoords.PSFIndexOfFacetCenter]
    
    plt.figure(584); plt.clf()
    fig, axes = plt.subplots(nrows=4, ncols=2, figsize=(12,10), num=584); 
    for ax,test in zip(axes.flatten(), fullSkyTests):
        s = allResults[test]['specifications']
        coords = allResults[test]['coords']
        #hp.mollview(np.log10(np.abs(allResults[test]['PSF'][coords.newPSFIndices,coords.PSFIndexOfFacetCenter])), title=test, fig=100+i, min=-5)
        plt.sca(ax)
        hp.mollview(np.log10(np.abs(allResults[test]['PSF'][coords.newPSFIndices,coords.PSFIndexOfFacetCenter])), title=test, min=-5, max = -2, cmap=pylab.cm.inferno, hold = True, xsize=500)
  

#%%
if False:
    HighResPSFTests = ['HERA-331 High-Res',
#                    'HERA-331 + Fidicuial Outriggers High-Res',
                    'HERA-331 + Redundant Outriggers High-Res',
                    'HERA-331 + Fidicuial Outriggers + Fidicuial Inriggers High-Res',
                    'HERA-331 + Redundant Outriggers + 6 Corner Inrigger Pairs High-Res',
                    'HERA-331 Split Core + Split Outriggers High-Res',
                    'HERA-331 Nine-Way Split Core + Split Outriggers High-Res']
    
    def logPSFComp(tests, figNum, nRows=2, nCols=3, minLog=-3, xSize=1600):
        plt.figure(figNum, figsize=(14,10)); plt.clf()
        fig, axes = plt.subplots(nrows=nRows, ncols=nCols, num=figNum, figsize=(14,10)); 
        for ax,test in zip(axes.flatten(), tests):
            s = allResults[test]['specifications']
            PSFIndices = hp.query_disc(s.mapNSIDE, hp.ang2vec(np.pi/2, 0), s.facetSize * 2*np.pi/360.0 * s.PSFextensionBeyondFacetFactor)            
            coords = allResults[test]['coords']
            plt.sca(ax)
            mapToPlot = np.zeros(coords.mapPixels)
            mapToPlot[PSFIndices] = allResults[test]['PSF'][coords.facetIndexOfFacetCenter,:]    
            hp.gnomview(np.log10(np.abs(mapToPlot)), title=test.replace('+','+\n'), max=0, min=minLog, hold=True, xsize=xSize)
    
    logPSFComp(HighResPSFTests, 7, xSize=200, minLog=-4)
    logPSFComp(HighResPSFTests, 8, xSize=800, minLog=-3)
    logPSFComp(HighResPSFTests, 9, xSize=1600, minLog=-2)

#def plotFacet(s,coords,facetMap,plotTitle):
#    if s.makeFacetSameAsAdaptivePSF and s.useAdaptiveHEALPixForPSF:
#        hp.mollview(facetMap[coords.newPSFIndices], title=plotTitle)
#    else:
#        mapToPlot = np.zeros(coords.mapPixels)
#        mapToPlot[coords.facetIndices] = facetMap    
#        hp.mollview(mapToPlot, title=plotTitle)
#        plt.axis(np.asarray([-1,1,-1,1]) * s.facetSize/50)



#%% Plotting
#
#plt.close('all')
#for i,folder in enumerate(ResultsToExamine):
#    print "Now plotting results from " + folder
#    s = allSpecs[i]
#    
#    plt.figure(figsize=(14,6))
#    plt.subplot(1, 2, 1)    
#
#    #plt.imshow(AtransA, interpolation='none')
#    #plt.imshow(np.linalg.pinv(AtransA)[0:s.nAntennas,0:s.nAntennas], interpolation='none')   
#  #  plt.scatter(s.antennaPositions[:,0], s.antennaPositions[:,1], c=gainVariances, s=100, norm=matplotlib.colors.LogNorm())
#    plt.scatter(s.antennaPositions[:,0], s.antennaPositions[:,1], c=allGainVariances[i], s=100)
#    plt.colorbar()
#    print "AtransA has " + str(allNZeroEVs[i][0]) + " zero eigenvalues."
#    print "BtransB has " + str(allNZeroEVs[i][1]) + " zero eigenvalues."
#    datacursor(display='single',formatter="x={x:.2f}\ny={y:.2f}".format)
#
#    plt.subplot(1, 2, 2)
#    plt.semilogy(np.real(allnoiseEvals[i]))
#    plt.xlabel('Eigenvalue Number',fontsize=16)
#    plt.ylabel('Eigenvalue of $\mathbf{A}^\dagger\mathbf{N}^{-1}\mathbf{A}$',fontsize=16)
#    #datacursor(display='single',formatter="N={x:.2f}\nEV={y:.2e}".format)
#
#
##%% Compare Hex Inriggers to non-hex
#
#allSkyHexIndex = ResultsDict['HERA331_and_Hex_inriggers_low_res_full_sky']
#allSkyInriggerIndex = ResultsDict['HERA331_and_3_inriggers_low_res_full_sky']
#plt.figure()
#plt.semilogy(allnoiseEvals[allSkyInriggerIndex])
#plt.semilogy(allnoiseEvals[allSkyHexIndex])
#plt.xlabel('Eigenvalue Number',fontsize=16)
#plt.ylabel('Eigenvalue of $\mathbf{A}^\dagger\mathbf{N}^{-1}\mathbf{A}$',fontsize=16)
#plt.legend(['3 Single Inriggers', '3 Hex Inriggers'])


#%% Array Figures

tests = ['HERA-331 + Fiducial Inriggers Full Sky',
         'HERA-331 + 6 Corner Inrigger Pairs Full Sky',
         'HERA-331 Split Core Full Sky']


plt.figure(101, figsize=(6,14)); plt.clf()

fig, axes = plt.subplots(nrows=3, ncols=1, num=101)
plt.subplots_adjust(hspace = 0.0)
for ax,test,i in zip(axes.flatten(), tests, range(3)):
    s = allResults[test]['specifications']
    hexScatter = dishScatter(fig, ax, s.antennaPositions[:,0], s.antennaPositions[:,1], np.ones(len(s.antennaPositions)), pylab.cm.gray, 7.0, cLim=[.6,2.0], blackOutlines = 1)
    ax.set_ylabel('Position (m)',size=16)
    if i < 2: ax.set_xticklabels([])
    else: ax.set_xlabel('Position (m)',size=16)
    
    ax.set_ylim([-180,180])        
    ax.set_xlim([-180,180])        
    ax.text(-150,150,['(a)','(b)','(c)'][i],size=20,ha='center',va='center')
    ax.set_aspect('equal', adjustable='box')
plt.tight_layout()
plt.subplots_adjust(hspace = 0.0)
fig.savefig('/Users/jsdillon/Desktop/Outrigger_Project/Paper/Figures/arrays.eps', format='eps')

#%% Simulaneous Redundancy

tests = ['HERA-331 + Fiducial Inriggers Full Sky',
         'HERA-331 + 6 Corner Inrigger Pairs Full Sky',
         'HERA-331 Split Core Full Sky',
         #'HERA-331 Nine-Way Split Core Full Sky'
         ]

plt.figure(102, figsize=(12,14)); plt.clf()

fig, axes = plt.subplots(nrows=len(tests), ncols=2, num=102, sharex = 'col')
plt.subplots_adjust(hspace = 0.0)
for i,ax in enumerate(axes.flatten()):
    s = allResults[tests[i/2]]['specifications']
    allBaselines = np.append(s.baselines,-s.baselines,axis=0)
    allRedundancies = np.append(s.baselineRedundancies,s.baselineRedundancies,axis=0)
    

    
    if i%2 == 0:
        scatterPlot = ax.scatter(allBaselines[:,0],allBaselines[:,1],c=allRedundancies,s=12,vmin=0, vmax=300,cmap=pylab.cm.inferno, edgecolors='none') 
        ax.set_ylim([-310,310])        
        ax.set_xlim([-310,310])        
        ax.text(-265,265,['(a)','(b)','(c)','(d)'][i/2],size=20,ha='center',va='center')
        ax.set_ylabel('Baseline (m)',size=16)    
    else: 
        ax.scatter(allBaselines[:,0],allBaselines[:,1],c=allRedundancies,s=60,vmin=0, vmax=300, cmap=pylab.cm.inferno, linewidths=.25, edgecolors='k')
        ax.set_ylim([-125,125])        
        ax.set_xlim([-125,125])
#       ax.set_yticklabels([])
    
    if i >= 4: ax.set_xlabel('Baseline (m)',size=16)
    
    ax.set_aspect('equal', adjustable='box-forced')
    
plt.tight_layout()
    
plt.subplots_adjust(hspace = 0.02, wspace = .1, right=0.85, left = .07)
cax = plt.axes([.87, axes.flat[5].get_position().y0, .05, axes.flat[0].get_position().y0+axes.flat[0].get_position().height - axes.flat[5].get_position().y0])
clr = fig.colorbar(scatterPlot, cax=cax, cmap = pylab.cm.jet.name)
clr.set_label('Simultaneous Redundancy', size=16)
 
fig.savefig('/Users/jsdillon/Desktop/Outrigger_Project/Paper/Figures/redundancy.eps', format='eps')


#%% Outriggers 

tests = ['HERA-331 + Fidicuial Outriggers + Fidicuial Inriggers High-Res',
         'HERA-331 + Redundant Outriggers + 6 Corner Inrigger Pairs High-Res',
         'HERA-331 Split Core + Split Outriggers High-Res']


plt.figure(103, figsize=(6,14)); plt.clf()

fig, axes = plt.subplots(nrows=3, ncols=1, num=103)
plt.subplots_adjust(hspace = 0.0)
for ax,test,i in zip(axes.flatten(), tests, range(3)):
    s = allResults[test]['specifications']
    hexScatter = dishScatter(fig, ax, s.antennaPositions[:,0], s.antennaPositions[:,1], np.ones(len(s.antennaPositions)), pylab.cm.gray, 7.0, cLim=[.6,2.0], blackOutlines = .25)
    ax.set_ylabel('Position (m)',size=16)
    if i < 2: ax.set_xticklabels([])
    else: ax.set_xlabel('Position (m)',size=16)
    
    ax.set_ylim([-550,550])        
    ax.set_xlim([-550,550])        
    ax.text(-450,450,['(a)','(b)','(c)'][i],size=20,ha='center',va='center')
    ax.set_aspect('equal', adjustable='box')
plt.tight_layout()
plt.subplots_adjust(hspace = 0.0)
fig.savefig('/Users/jsdillon/Desktop/Outrigger_Project/Paper/Figures/arraysFull.eps', format='eps')


#%% Simulaneous Redundancy Full

tests = ['HERA-331 + Fidicuial Outriggers + Fidicuial Inriggers High-Res',
         'HERA-331 + Redundant Outriggers + 6 Corner Inrigger Pairs High-Res',
         'HERA-331 Split Core + Split Outriggers High-Res']


plt.figure(104, figsize=(12,14)); plt.clf()

fig, axes = plt.subplots(nrows=3, ncols=2, num=104, sharex = 'col')
plt.subplots_adjust(hspace = 0.0)
for i,ax in enumerate(axes.flatten()):
    s = allResults[tests[i/2]]['specifications']
    allBaselines = np.append(s.baselines,-s.baselines,axis=0)
    allRedundancies = np.append(s.baselineRedundancies,s.baselineRedundancies,axis=0)
    

    
    if i%2 == 0:
        scatterPlot = ax.scatter(allBaselines[:,0],allBaselines[:,1],c=allRedundancies,edgecolors='none',s=5,vmin=0, vmax=300, cmap=pylab.cm.inferno, linewidths=.5)
        #plt.set_clim([0,300])    
        ax.set_ylim([-550,550])        
        ax.set_xlim([-550,550])        
        ax.text(-475,475,['(a)','(b)','(c)'][i/2],size=20,ha='center',va='center')
        ax.set_ylabel('Baseline (m)',size=16)    
    else: 
        ax.scatter(allBaselines[:,0],allBaselines[:,1],c=allRedundancies,s=30,vmin=0, vmax=300, cmap=pylab.cm.inferno, linewidths=.25, edgecolors='k')        
        ax.set_ylim([0,399])        
        ax.set_xlim([0,399])
#       ax.set_yticklabels([])
    
    if i >= 4: ax.set_xlabel('Baseline (m)',size=16)
    
    ax.set_aspect('equal', adjustable='box-forced')
    
plt.tight_layout()
    
plt.subplots_adjust(hspace = 0.02, wspace = .1, right=0.85, left = .07)
cax = plt.axes([.87, axes.flat[5].get_position().y0, .05, axes.flat[0].get_position().y0+axes.flat[0].get_position().height - axes.flat[5].get_position().y0])
clr = fig.colorbar(scatterPlot, cax=cax, cmap = pylab.cm.jet.name)
clr.set_label('Simultaneous Redundancy', size=16)
 
fig.savefig('/Users/jsdillon/Desktop/Outrigger_Project/Paper/Figures/redundancyFull.eps', format='eps')

#%% Redundant Calibratability
test = 'HERA-331 Full Sky'

fig = plt.figure(105, figsize=(6.5,5)); plt.clf();
ax = plt.gca()
s = allResults[test]['specifications']
dishScatter(fig, ax, s.antennaPositions[:,0], s.antennaPositions[:,1], np.sqrt(allResults[test]['gainVariances']), pylab.cm.inferno, 7.0, blackOutlines=.5)
ax.set_ylim([-180,180])        
ax.set_xlim([-180,180])        
ax.set_xlabel('Position (m)',size=14)
ax.set_ylabel('Position (m)',size=14)
divider = make_axes_locatable(ax)    
cax = divider.append_axes("right", size="5%", pad=0.05)    
clr = plt.colorbar(cax=cax, format='%.4f')
clr.set_label('Relative Gain Calibration Error', size=14)
ax.set_aspect('equal', adjustable='box-forced')

plt.tight_layout()

fig.savefig('/Users/jsdillon/Desktop/Outrigger_Project/Paper/Figures/redCalCore.eps', format='eps')

#plt.figure(1000)
#plt.plot(gainVariances)

#%% Redundant Calibratability (all arrays)
tests = [['HERA-331 + Redundant Outriggers + 6 Corner Inrigger Pairs High-Res',
          'HERA-331 Split Core + Split Outriggers High-Res'],
          ['HERA-331 + 6 Corner Inrigger Pairs Full Sky',
          'HERA-331 Split Core Full Sky']]


for zoom in range(2):      
    plt.figure(106+zoom, figsize=(5.5,7)); plt.clf()
    fig, axes = plt.subplots(nrows=2, ncols=1, num=106+zoom, sharex = 'col')
    plt.subplots_adjust(hspace = 0.02, right=.7, left=.15,top=.97)
    for ax,test,i in zip(axes.flatten(), tests[zoom], range(6)):
        s = allResults[test]['specifications']
        thisScatter = dishScatter(fig, ax, s.antennaPositions[:,0], s.antennaPositions[:,1], 
                                  np.sqrt(allResults[test]['gainVariances']), pylab.cm.inferno, 7.0, blackOutlines = .25+.25*zoom, cLim = [(.0558,.1472),(.0556,.1139)][zoom])
        print test, np.min(np.sqrt(allResults[test]['gainVariances'])), np.max(np.sqrt(allResults[test]['gainVariances']))
        if zoom:
            ax.set_ylim([-180,180])        
            ax.set_xlim([-180,180])          
        else:
            ax.set_ylim([-550,550])        
            ax.set_xlim([-550,550])              
        if i == 1: ax.set_xlabel('Position (m)',size=14)
        ax.set_ylabel('Position (m)',size=14)
    
        if zoom: ax.text(-150,150,['(b)','(c)'][i],size=20,ha='center',va='center')
        else: ax.text(-450,450,['(b)','(c)'][i],size=20,ha='center',va='center')
    
        divider = make_axes_locatable(ax)    
#        cax = divider.append_axes("right", size="5%", pad=0.05)    
        #clr = plt.colorbar(cax=cax, format='%.4f')
        #clr.set_label('Relative Gain Calibration Error', size=14)
        ax.set_aspect('equal', adjustable='box-forced')

    cax = plt.axes([.75, axes.flat[1].get_position().y0, .06, axes.flat[0].get_position().y0+axes.flat[0].get_position().height - axes.flat[1].get_position().y0])
    clr = fig.colorbar(thisScatter, cax=cax, cmap = pylab.cm.inferno,format='%.4f')#, ticks=cticks)
    clr.set_label('Relative Gain Calibration Error', size=14)
    
    
#    plt.tight_layout()
#    plt.subplots_adjust(wspace = .3)
    if zoom: fig.savefig('/Users/jsdillon/Desktop/Outrigger_Project/Paper/Figures/redCalArraysZoom.eps', format='eps')
    else: fig.savefig('/Users/jsdillon/Desktop/Outrigger_Project/Paper/Figures/redCalArrays.eps', format='eps')

#%% Full sky mapmaking

fullSkyTests = ['HERA-331 Full Sky',
                'HERA-331 Full Sky Long Integration',                
                'HERA-331 + Fiducial Inriggers Full Sky',
                'HERA-331 + Fiducial Inriggers Full Sky Long Integration',                
                'HERA-331 + 6 Corner Inrigger Pairs Full Sky',
                'HERA-331 + 6 Corner Inrigger Pairs Full Sky Long Integration',
                'HERA-331 Split Core Full Sky',
                'HERA-331 Split Core Full Sky Long Integration']
legendtext =['HERA-331 Core Only','...with Rotation Synthesis', 'Fiducial Core Configuration (a)','...with Rotation Synthesis', 'Redundant Corners Configuration Core (b)','...with Rotation Synthesis', 'Split-Core Configuration (c)','...with Rotation Synthesis']
fig = plt.figure(582, figsize=(12,5)); plt.clf()
for i,test in enumerate(fullSkyTests):
    s = allResults[test]['specifications']
    spectrum = np.real(allResults[test]['noiseEvals'][allResults[test]['noiseEvals'] > 0])
    thisPlot = plt.semilogy(np.arange(1,len(spectrum)+1),spectrum/ s.nAntennas**2 * 331**2, ls = ['-','--'][i%2], c=pylab.cm.inferno((i/2)*.25))

maxEVs = np.max(np.asarray([len(allResults[test]['noiseEvals']) for test in fullSkyTests]))
plt.xlabel('Eigenvalue Number',fontsize=16)
plt.ylabel('Eigenvalue of $\mathbf{A}^\dagger\mathbf{N}^{-1}\mathbf{A}$',fontsize=16)
plt.legend(legendtext, loc='upper right', fontsize=12)
plt.ylim([1e-3, 1e6])
plt.xlim([1, maxEVs])
plt.tight_layout()
fig.savefig('/Users/jsdillon/Desktop/Outrigger_Project/Paper/Figures/fullSkyEVs.eps', format='eps')


#%% High-res mapmaking

highResTests = ['HERA-331 High-Res',
                'HERA-331 High-Res Long Integration',
                'HERA-331 + Fidicuial Outriggers + Fidicuial Inriggers High-Res',
                'HERA-331 + Fidicuial Outriggers + Fidicuial Inriggers High-Res Long Integration',
                'HERA-331 + Redundant Outriggers + 6 Corner Inrigger Pairs High-Res',
                'HERA-331 + Redundant Outriggers + 6 Corner Inrigger Pairs High-Res Long Integration',
                'HERA-331 Split Core + Split Outriggers High-Res',
                'HERA-331 Split Core + Split Outriggers High-Res Long Integration']
legendtext =['HERA-331 Core Only','...with Rotation Synthesis', 'Fiducial Configuration with Outriggers (a)','...with Rotation Synthesis', 'Redundant Corners Configuration with Outriggers (b)','...with Rotation Synthesis', 'Split-Core Configuration with Outriggers (c)','...with Rotation Synthesis']
fig = plt.figure(583, figsize=(12,5)); plt.clf()
for i,test in enumerate(highResTests):
    s = allResults[test]['specifications']
    spectrum = np.real(allResults[test]['noiseEvals'][allResults[test]['noiseEvals'] > 0])
    thisPlot = plt.semilogy(np.arange(1,len(spectrum)+1),spectrum/ s.nAntennas**2 * 331**2, ls = ['-','--'][i%2], c=pylab.cm.inferno((i/2)*.25))

maxEVs = np.max(np.asarray([len(allResults[test]['noiseEvals']) for test in highResTests]))
plt.xlabel('Eigenvalue Number',fontsize=16)
plt.ylabel('Eigenvalue of $\mathbf{A}^\dagger\mathbf{N}^{-1}\mathbf{A}$',fontsize=16)
plt.legend(legendtext, loc='lower left', fontsize=12)
plt.ylim([1e-6, 1e4])
plt.xlim([1, maxEVs])
plt.tight_layout()
fig.savefig('/Users/jsdillon/Desktop/Outrigger_Project/Paper/Figures/highResEVs.eps', format='eps')


  
    
#%% High Res PSFs

highResTests = [#'HERA-331 High-Res',
                #'HERA-331 High-Res Long Integration',
                'HERA-331 + Fidicuial Outriggers + Fidicuial Inriggers High-Res',
                'HERA-331 + Fidicuial Outriggers + Fidicuial Inriggers High-Res Long Integration',
                'HERA-331 + Redundant Outriggers + 6 Corner Inrigger Pairs High-Res',
                'HERA-331 + Redundant Outriggers + 6 Corner Inrigger Pairs High-Res Long Integration',
                'HERA-331 Split Core + Split Outriggers High-Res',
                'HERA-331 Split Core + Split Outriggers High-Res Long Integration']

def PSFplot(figNum,filename,size,logmin,logmax, cticks,plotOffCenter=None, allTests = highResTests, cbarlabel='$log_{10}$|PSF Amplitude|', save=True):
    plt.figure(figNum, figsize=(5,8)); plt.clf()
    fig, axes = plt.subplots(nrows=3, ncols=2, num=figNum, figsize=(5,8), sharex = 'col'); 
    plt.subplots_adjust(hspace = 0.0, wspace = .0, bottom = .16, top=.98, left = .15)
    for i,ax,test in zip(range(6),axes.flatten(), allTests):
        s = allResults[test]['specifications']        
        coords = allResults[test]['coords']
        plt.sca(ax)
        if plotOffCenter is not None:
            mapToPlot = allResults[test]['PSF'][coords.newPSFIndices,plotOffCenter]            
        else:
            PSFIndices = hp.query_disc(s.mapNSIDE, hp.ang2vec(np.pi/2, 0), s.facetSize * 2*np.pi/360.0 * s.PSFextensionBeyondFacetFactor)            
            mapToPlot = np.zeros(coords.mapPixels)
            mapToPlot[PSFIndices] = allResults[test]['PSF'][coords.facetIndexOfFacetCenter,:]
        colorplot = ax.imshow(np.asarray([[logmin,logmax]]),cmap = pylab.cm.inferno)
        hp.cartview(np.log10(np.abs(mapToPlot)), title='', max=logmax, min=logmin, hold=True, cmap=pylab.cm.inferno, cbar=False, lonra=[-size,size], latra=[-size,size], xsize=2000)
        ax.set_aspect('equal', adjustable='box-forced')
    
    
        ax1 = fig.add_subplot(3,2,i+1)
        ax1.patch.set_alpha(0)
        ax1.set_ylim([-size,size])        
        ax1.set_xlim([-size,size])        
        ax1.set_aspect('equal', adjustable='box-forced')
        if i%2==0: 
            ax1.set_ylabel('Sky Position (deg)')
            ax1.text(-9.5,9.5,['(a)','(b)','(c)'][i/2],size=16,ha='center',va='center',color='w')
        else: ax1.set_yticklabels([])
        if i >= 4: ax1.set_xlabel('Sky Position (deg)')
        else: ax1.set_xticklabels([])
    
    cax = plt.axes([axes.flat[4].get_position().x0, .07, axes.flat[5].get_position().x0 + axes.flat[5].get_position().width - axes.flat[4].get_position().x0, .04])
    clr = fig.colorbar(colorplot, cax=cax, cmap = pylab.cm.inferno, orientation='horizontal', ticks=cticks)
    clr.set_label(cbarlabel, size=16)
    if save: fig.savefig('/Users/jsdillon/Desktop/Outrigger_Project/Paper/Figures/'+filename, format='eps', dpi=600)

PSFplot(585,'highResPSFs.eps',12.5,-4,0, [-4,-3.5,-3.0,-2.5,-2.0, -1.5, -1,-.5,0])



#%%
#NOTE TO SELF: Now I need to find representative pixels for this
#if False:
#    plt.figure(2323); plt.clf()
#    s = allResults[test]['specifications']        
#    coords = allResults[test]['coords']
#    mapToPlot = allResults[test]['PSF'][coords.newPSFIndices,coords.PSFIndexOfFacetCenter]
#    hp.mollview(np.log10(np.abs(mapToPlot)), fig=2323)
#    
#    plt.figure(2324); plt.clf()
#    mapToPlot = np.arange(len(allResults[test]['PSF']))[coords.newPSFIndices] -5540+24
#    hp.mollview(mapToPlot, fig=2324)

#%%


#
#fullSkyTests = ['HERA-331 Full Sky',
##                'HERA-331 + Fiducial Inriggers Full Sky',
##                'HERA-331 + 6 Corner Inrigger Pairs Full Sky',
#                'HERA-331 Split Core Full Sky']
#psPixels = [4530, 5564, 4587, 4676, 5380]
#
#
##PSFplot(685,'offCenter.eps',12.5,-5,0, [-5,-4,-3.0,-2.0, -1,0], plotOffCenter=1000, allTests = fullSkyTests, cbarlabel='$log_{10}$|Relative Flux|', save=False)
#s = allResults['HERA-331 Full Sky']['specifications']
#coords = allResults['HERA-331 Full Sky']['coords']
#Geometry.convertEquatorialToHorizontal(s,coords.PSFRAs[0],coords.PSFDecs[0],s.facetRA*24.0/360.0)
#for pix in range(coords.nPSFPixels):
#    alt =  Geometry.convertEquatorialToHorizontal(s,coords.PSFRAs[pix],coords.PSFDecs[pix],s.facetRA*24.0/360.0)[0]*360.0/2/np.pi
#    if alt < 1.5 and alt > .5: print pix, alt
#
#
#logmin = -2
#logmax = 2
#size = 12.5
#
#plt.figure(685, figsize=(12,14)); plt.clf()
#fig, axes = plt.subplots(nrows=len(psPixels), ncols=len(fullSkyTests), num=685, figsize=(12,14), sharex = 'col'); 
#plt.subplots_adjust(hspace = 0.0, wspace = .0, bottom = .16, top=.98, left = .15)
#
#    
#for i,test in enumerate(fullSkyTests):
#    for j,psPix in enumerate(psPixels):
#        ax = axes[j][i]
#        s = allResults[test]['specifications']        
#        coords = allResults[test]['coords']
#        plt.sca(ax)
#    
##        mapToPlot = allResults[test]['PSF'][psPix,coords.newPSFIndices]# - allResults[fullSkyTests[0]]['PSF'][psPix,coords.newPSFIndices]
#        thisPSF = allResults[test]['PSF'][psPix,coords.newPSFIndices]
#        firstPSF = allResults[fullSkyTests[0]]['PSF'][psPix,coords.newPSFIndices]
#        mapToPlot = thisPSF / firstPSF
#        colorplot = ax.imshow(np.asarray([[logmin,logmax]]),cmap = pylab.cm.inferno)
#        hp.cartview(np.log10(np.abs(mapToPlot)), title='', max=logmax, min=logmin, hold=True, cmap=pylab.cm.inferno, cbar=False, lonra=[-size,size], latra=[-size,size])#, xsize=2000)
#        print np.max(np.log10(np.abs(mapToPlot)))        
#        ax.set_aspect('equal', adjustable='box-forced')
#    #
#    #
#        ax1 = fig.add_subplot(len(psPixels),len(fullSkyTests),i*len(psPixels)+j+1)
#        ax1.patch.set_alpha(0)
#        ax1.set_ylim([-size,size])        
#        ax1.set_xlim([-size,size])        
#        ax1.set_aspect('equal', adjustable='box-forced')
#    #    if i%2==0: 
#    #        ax1.set_ylabel('Sky Position (deg)')
#    #        ax1.text(-9.5,9.5,['(a)','(b)','(c)'][i/2],size=16,ha='center',va='center',color='w')
#    #    else: ax1.set_yticklabels([])
#    #    if i >= 4: ax1.set_xlabel('Sky Position (deg)')
#    #    else: ax1.set_xticklabels([])
#
##cax = plt.axes([axes.flat[4].get_position().x0, .07, axes.flat[5].get_position().x0 + axes.flat[5].get_position().width - axes.flat[4].get_position().x0, .04])
##clr = fig.colorbar(colorplot, cax=cax, cmap = pylab.cm.inferno, orientation='horizontal', ticks=cticks)
##clr.set_label(cbarlabel, size=16)
##if save: fig.savefig('/Users/jsdillon/Desktop/Outrigger_Project/Paper/Figures/'+filename, format='eps', dpi=600)


#%% PSF Inversion

#if False:
#    tests = ['HERA-331 Full Sky Long Integration',                
#             'HERA-331 + Fiducial Inriggers Full Sky Long Integration',                
#             'HERA-331 + 6 Corner Inrigger Pairs Full Sky Long Integration',
#             'HERA-331 Split Core Full Sky Long Integration']
#    #tests = ['HERA-331 Full Sky']#,                
#    #         'HERA-331 + Fiducial Inriggers Full Sky',                
#    #         'HERA-331 + 6 Corner Inrigger Pairs Full Sky',
#    #         'HERA-331 Split Core Full Sky']
#    Pinvs = {}
#    for i,test in enumerate(tests):
#        print "Working on " + test
#        thisTest = allResults[test]            
#        D = thisTest['Dmatrix']
#    #    Pinvs[test] = np.linalg.pinv(np.diag(D**-1).dot(PSF), rcond = 1e-12)
#        Pinvs[test] = np.diag(D).dot(PSF)
#    
#    #%%
#    for i,test in enumerate(tests):    
#        plt.figure(2323+i); plt.clf()
#        thisTest = allResults[test]            
#        s = thisTest['specifications']
#        coords = thisTest['coords']
#        mapToPlot = (np.sqrt(np.diag(Pinvs[test]))*(coords.newPSFNSIDEs/np.max(coords.newPSFNSIDEs))**2)[coords.newPSFIndices] 
#        mapToPlot2 = (np.sqrt(np.diag(Pinvs[tests[0]]))*(coords.newPSFNSIDEs/np.max(coords.newPSFNSIDEs))**2)[coords.newPSFIndices] 
#        size = 89    
#        hp.cartview(np.abs(mapToPlot), fig=2323+i,lonra=[-size,size], latra=[-size,size], cmap=pylab.cm.inferno)


#%% Full Sky Eigenvector Plots
if True:
    
                        #[np.arange(0,100), 
                        #np.arange(0,630), 
                        #np.arange(0,1000),
                        #np.arange(0,1200),
                        #np.arange(1200,5900),
                        #np.arange(1300,1500),
                        #np.arange(630,3000),
                        #np.arange(1000,2000), 
                        #np.arange(1000,3000),
                        #np.arange(1000,4000),
                        #np.arange(3000,4000),
                        #np.arange(4000,5900),
                        #np.arange(4150,4750),
                        #np.arange(4750,5900)]
    rangesOfInterest = [[(0,629)],[(630,2999),(3000,3999),(4000,5899)]]
    fullSkyTests = [#'HERA-331 Full Sky',
                    #'HERA-331 Full Sky Long Integration',                
                    'HERA-331 + Fiducial Inriggers Full Sky',
                    #'HERA-331 + Fiducial Inriggers Full Sky Long Integration',                
                    'HERA-331 + 6 Corner Inrigger Pairs Full Sky',
                    #'HERA-331 + 6 Corner Inrigger Pairs Full Sky Long Integration',
                    'HERA-331 Split Core Full Sky']#,
                    #'HERA-331 Split Core Full Sky Long Integration']
    logmin,logmax = -3,3
    size = 89
    
    plt.figure(172, figsize=(11,14)); plt.clf()
    fig, axes = plt.subplots(ncols=len(rangesOfInterest), nrows=len(fullSkyTests), num=172, figsize=(12,14), sharex = 'col'); 
    plt.subplots_adjust(hspace = 0.02, wspace = .02, left = .07, right = .84, top = .94, bottom = .06)

    for i,test in enumerate(fullSkyTests):
        result = allResults[test]
        coords = result['coords']
        Dmatrix = result['Dmatrix']
        allResults[test]['specifications']        
        PSF = result['PSF']
        evecLen = len(result['noiseEvecs'][0])
        posVarIndicies = np.where(np.diag(PSF[:,coords.facetIndexLocationsInPSFIndexList])>0)        
        for j,rng in enumerate(rangesOfInterest):
            evalsToPlot = np.zeros(len(coords.PSFDecs))
            for subrng in rng:            
                evalsToPlot[posVarIndicies] += result['partialEigenspace_'+str(subrng[0])+'-'+str(subrng[-1])] 
            evalsToPlot *= (coords.newPSFNSIDEs/np.max(coords.newPSFNSIDEs))**4 / s.nAntennas**2 * 331**2 
            ax = axes[i][j]
            s = allResults[test]['specifications']        
            plt.sca(ax)
            colorplot = ax.imshow(np.asarray([[logmin,logmax]]),cmap = pylab.cm.inferno)
            hp.cartview(np.log10(np.abs(evalsToPlot[coords.newPSFIndices])), lonra=[-size,size], latra=[-size,size], max=logmax, min=logmin, hold=True, cmap=pylab.cm.inferno,cbar=False,title='')
   
            ax.set_aspect('equal', adjustable='box-forced')
            ax1 = fig.add_subplot(len(fullSkyTests),len(rangesOfInterest),i*len(rangesOfInterest)+j+1)
            ax1.patch.set_alpha(0)
            ax1.set_ylim([-size,size])        
            ax1.set_xlim([-size,size])        
            ax1.set_aspect('equal', adjustable='box-forced')
            ax1.text(-size*.85,size*.85,['(a)','(b)','(c)'][i],size=18,ha='center',va='center',color='w')
            thisTitle = ''
            if i == 0: thisTitle = ['$N_\lambda \leq 630$','$N_\lambda > 630$'][j]
            ax1.set_title(thisTitle, size=18)
            if j==0: 
                ax1.set_ylabel('Sky Position (deg)',size=16)
                #ax1.text(-9.5,9.5,['(a)','(b)','(c)'][i/2],size=16,ha='center',va='center',color='w')
            else: ax1.set_yticklabels([])
            if i >= 2: ax1.set_xlabel('Sky Position (deg)',size=16)
            else: ax1.set_xticklabels([])
            [t.set_color('w') for t in ax1.xaxis.get_ticklines()]
            for child in ax1.get_children():
                if isinstance(child, matplotlib.spines.Spine):
                    child.set_color('w')
    
    cax = plt.axes([.87, axes.flat[5].get_position().y0, .04, axes.flat[0].get_position().y0+axes.flat[0].get_position().height - axes.flat[5].get_position().y0])

    clr = fig.colorbar(colorplot, cax=cax, cmap = pylab.cm.inferno)#, ticks=cticks)
    clr.set_label('$log_{10}[$diagonal of the projected $\mathbf{A}^\dagger\mathbf{N}^{-1}\mathbf{A}]$', size=16)
    if True: fig.savefig('/Users/jsdillon/Desktop/Outrigger_Project/Paper/Figures/fullSkyPartialEigenspaces.eps', format='eps', dpi=600)


#%% High Res Eigenvector Plots
if True:
    
                        #[np.arange(0,100), 
                        #np.arange(0,630), 
                        #np.arange(0,1000),
                        #np.arange(0,1200),
                        #np.arange(1200,5900),
                        #np.arange(1300,1500),
                        #np.arange(630,3000),
                        #np.arange(1000,2000), 
                        #np.arange(1000,3000),
                        #np.arange(1000,4000),
                        #np.arange(3000,4000),
                        #np.arange(4000,5900),
                        #np.arange(4150,4750),
                        #np.arange(4750,5900)]
    rangesOfInterest = [[(0,999)],[(1000,3999),(3000,3999),(4000,5899)]]
    highResTests = [#'HERA-331 High-Res',
                    #'HERA-331 High-Res Long Integration',
                    'HERA-331 + Fidicuial Outriggers + Fidicuial Inriggers High-Res',
  #                  'HERA-331 + Fidicuial Outriggers + Fidicuial Inriggers High-Res Long Integration',
                    'HERA-331 + Redundant Outriggers + 6 Corner Inrigger Pairs High-Res',
 #                   'HERA-331 + Redundant Outriggers + 6 Corner Inrigger Pairs High-Res Long Integration',
                    'HERA-331 Split Core + Split Outriggers High-Res']#,
#                    'HERA-331 Split Core + Split Outriggers High-Res Long Integration']

    
    logmin,logmax = -1,1
    size = 10
    
    plt.figure(173, figsize=(11,14)); plt.clf()
    fig, axes = plt.subplots(ncols=len(rangesOfInterest), nrows=len(highResTests), num=173, figsize=(12,14), sharex = 'col'); 
    plt.subplots_adjust(hspace = 0.02, wspace = .02, left = .07, right = .84, top = .94, bottom = .06)

    for i,test in enumerate(highResTests):
        result = allResults[test]
        coords = result['coords']
        Dmatrix = result['Dmatrix']
        allResults[test]['specifications']        
        PSF = result['PSF']
        evecLen = len(result['noiseEvecs'][0])
        for j,rng in enumerate(rangesOfInterest):
            evalsToPlot = np.zeros(len(coords.facetDecs))
            for subrng in rng:            
                evalsToPlot += result['partialEigenspace_'+str(subrng[0])+'-'+str(subrng[-1])] 
            evalsToPlot *= 331.0**2/s.nAntennas**2
            ax = axes[i][j]
            s = allResults[test]['specifications']        
            plt.sca(ax)
            colorplot = ax.imshow(np.asarray([[logmin,logmax]]),cmap = pylab.cm.inferno)
            mapToPlot = np.zeros(coords.mapPixels)
            mapToPlot[coords.facetIndices] = evalsToPlot
            
            
            hp.cartview(np.log10(np.abs(mapToPlot)),lonra=[-size,size], latra=[-size,size], max=logmax, min=logmin, hold=True, cmap=pylab.cm.inferno,cbar=False,title='')
   
            ax.set_aspect('equal', adjustable='box-forced')
            ax1 = fig.add_subplot(len(highResTests),len(rangesOfInterest),i*len(rangesOfInterest)+j+1)
            ax1.patch.set_alpha(0)
            ax1.set_ylim([-size,size])        
            ax1.set_xlim([-size,size])        
            ax1.set_aspect('equal', adjustable='box-forced')
            ax1.text(-size*.85,size*.85,['(a)','(b)','(c)'][i],size=18,ha='center',va='center',color='w')
            thisTitle = ''
            if i == 0: thisTitle = ['$N_\lambda \leq 1000$','$N_\lambda > 1000$'][j]
            ax1.set_title(thisTitle, size=18)
            if j==0: 
                ax1.set_ylabel('Sky Position (deg)',size=16)
                #ax1.text(-9.5,9.5,['(a)','(b)','(c)'][i/2],size=16,ha='center',va='center',color='w')
            else: ax1.set_yticklabels([])
            if i >= 2: ax1.set_xlabel('Sky Position (deg)',size=16)
            else: ax1.set_xticklabels([])
            [t.set_color('w') for t in ax1.xaxis.get_ticklines()]
            for child in ax1.get_children():
                if isinstance(child, matplotlib.spines.Spine):
                    child.set_color('w')
    
    cax = plt.axes([.87, axes.flat[5].get_position().y0, .04, axes.flat[0].get_position().y0+axes.flat[0].get_position().height - axes.flat[5].get_position().y0])

    clr = fig.colorbar(colorplot, cax=cax, cmap = pylab.cm.inferno)#, ticks=cticks)
    clr.set_label('$log_{10}[$diagonal of the projected $\mathbf{A}^\dagger\mathbf{N}^{-1}\mathbf{A}]$', size=16)
    if True: fig.savefig('/Users/jsdillon/Desktop/Outrigger_Project/Paper/Figures/highResPartialEigenspaces.eps', format='eps', dpi=600)
