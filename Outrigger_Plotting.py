import numpy as np
import healpy as hp
import math
import sys,os,os.path
import matplotlib
import matplotlib.pyplot as plt
from MAPS21cm.Specifications import Specifications
from MAPS21cm import Geometry
import cPickle as pickle
import scipy.constants as const
from mpldatacursor import datacursor
from os import listdir
from os.path import isfile, join
import pylab
plt.close('all')

scriptDirectory = os.path.dirname(os.path.abspath(__file__))
allTests = pickle.load(open(scriptDirectory + "/Results/allTests.p","rb"))

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


#%% Load results of calcuations
allResults = {}
for test in allTests:
    allResults[test['title']] = {}
    resultsDir =  scriptDirectory + '/Results/' + test['folder']
    onlyFiles = [ f for f in listdir(resultsDir) if isfile(join(resultsDir,f)) ]
    print 'Now loading ' + test['title']
    for f in onlyFiles:      
        if f.split('.')[1] == 'npy': allResults[test['title']][f.split('.')[0]] = np.load(resultsDir + f)
        if f.split('.')[1] == 'p': allResults[test['title']][f.split('.')[0]] = pickle.load(open(resultsDir + f,'rb'))
        if f.split('.')[1] == 'txt': allResults[test['title']][f.split('.')[0]] = np.loadtxt(resultsDir + f)

#%% Redundant Calibratability
test = 'HERA-331 Full Sky'

def dishScatter(fig, ax, xPos, yPos, cVals, colormap, radii=7.0, cLim = None): 
    if not hasattr(radii,'len'): radii = np.ones(len(xPos))*radii    
    if cLim is None:     
        underlyingScatter = plt.scatter(xPos, yPos, c=cVals, s=0, cmap=colormap.name)
        cVals = (cVals - np.min(cVals))/(np.max(cVals) - np.min(cVals))
    else:
        underlyingScatter = plt.scatter(xPos, yPos, c=cVals, s=0, cmap=colormap.name, vmin = cLim[0], vmax = cLim[1])
        cVals = 1.0*(cVals - cLim[0])/(cLim[1]-cLim[0])
    for x,y,r,c in zip(xPos, yPos, radii, colormap(cVals)): 
        ax.add_patch(plt.Circle((x,y), r, fc=c, ec=c)) 
    return underlyingScatter

fig = plt.figure(1); plt.clf();
ax = plt.gca()
s = allResults[test]['specifications']
dishScatter(fig, ax, s.antennaPositions[:,0], s.antennaPositions[:,1], np.sqrt(allResults[test]['gainVariances']), pylab.cm.jet, 7.0)
plt.colorbar(format='%.3f')
plt.axis('equal')    
plt.tight_layout()


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
         'HERA-331 + 3 Corner Inrigger Pairs Full Sky',
         'HERA-331 + 6 Corner Inrigger Pairs Full Sky',
         'HERA-331 Split Core Full Sky']
CalibratabilityComparison(inriggerTests,2)


fullArrayTests = ['HERA-331 + Fidicuial Outriggers + Fidicuial Inriggers High-Res',
         'HERA-331 + Redundant Outriggers + 3 Corner Inrigger Pairs High-Res',
         'HERA-331 + Redundant Outriggers + 6 Corner Inrigger Pairs High-Res',
         'HERA-331 Split Core + Split Outriggers High-Res']
CalibratabilityComparison(fullArrayTests,3)      


#%% Full sky mapmaking

fullSkyTests = ['HERA-331 Full Sky',
                'HERA-331 + Fiducial Inriggers Full Sky',
                'HERA-331 + 3 Corner Inrigger Pairs Full Sky',
                'HERA-331 + 6 Corner Inrigger Pairs Full Sky',
                'HERA-331 Split Core Full Sky']
plt.figure(4, figsize=(12,8)); plt.clf()
for test in fullSkyTests:
    plt.semilogy(np.real(allResults[test]['noiseEvals']))
plt.xlabel('Eigenvalue Number',fontsize=16)
plt.ylabel('Eigenvalue of $\mathbf{A}^\dagger\mathbf{N}^{-1}\mathbf{A}$',fontsize=16)
plt.legend(fullSkyTests, loc='upper right', fontsize=11)
plt.ylim([1e-4, 1e6])
plt.title('Full Sky Mapmaking')

#%% High res mapmaking

highResTests = ['HERA-331 High-Res',
                'HERA-331 + Fidicuial Outriggers High-Res',
                'HERA-331 + Redundant Outriggers High-Res',
                'HERA-331 + Fidicuial Outriggers + Fidicuial Inriggers High-Res',
                'HERA-331 + Redundant Outriggers + 3 Corner Inrigger Pairs High-Res',
                'HERA-331 + Redundant Outriggers + 6 Corner Inrigger Pairs High-Res',
                'HERA-331 Split Core + Split Outriggers High-Res']
plt.figure(5, figsize=(12,8)); plt.clf()
for test in highResTests:
    plt.semilogy(np.real(allResults[test]['noiseEvals']))
plt.xlabel('Eigenvalue Number',fontsize=16)
plt.ylabel('Eigenvalue of $\mathbf{A}^\dagger\mathbf{N}^{-1}\mathbf{A}$',fontsize=16)
plt.legend(highResTests, loc='upper right', fontsize=11)
plt.ylim([1e-6, 1e4])
plt.title('High Res Mapmaking')

#%% Full sky PSFs

fullSkyCoords = allResults['HERA-331 Full Sky']['coords']
HERAFullSkyPSF = allResults['HERA-331 Full Sky']['PSF'][fullSkyCoords.newPSFIndices,fullSkyCoords.PSFIndexOfFacetCenter]

plt.figure(6); plt.clf()
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12,10), num=6); 
for ax,test in zip(axes.flatten(), fullSkyTests[1:]):
    s = allResults[test]['specifications']
    coords = allResults[test]['coords']
    #hp.mollview(np.log10(np.abs(allResults[test]['PSF'][coords.newPSFIndices,coords.PSFIndexOfFacetCenter])), title=test, fig=100+i, min=-5)
    plt.sca(ax)
    hp.mollview(np.log10(np.abs(HERAFullSkyPSF)) - np.log10(np.abs(allResults[test]['PSF'][coords.newPSFIndices,coords.PSFIndexOfFacetCenter])), 
                title=('Difference of Log10 of \n' + test), min=-2, max = 2, cmap='RdBu', hold=True)


#%%
HighResPSFTests = ['HERA-331 High-Res',
                'HERA-331 + Fidicuial Outriggers High-Res',
                'HERA-331 + Redundant Outriggers High-Res',
                'HERA-331 + Fidicuial Outriggers + Fidicuial Inriggers High-Res',
                'HERA-331 + Redundant Outriggers + 6 Corner Inrigger Pairs High-Res',
                'HERA-331 Split Core + Split Outriggers High-Res']

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
