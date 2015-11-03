import numpy as np
import healpy as hp
import math
import sys,os,os.path
import matplotlib
import matplotlib.pyplot as plt
from MAPS21cm.Specifications import Specifications
from MAPS21cm import Geometry
import scipy.constants as const
from mpldatacursor import datacursor

ResultsToExamine = ['HERA331',
                    'HERA331_and_3_inriggers',
                    'HERA331_and_Joshs_outriggers',
                    'HERA331_and_Joshs_outriggers_and_3_inriggers',
                    'HERA331_low_res_full_sky',
                    'HERA331_and_3_inriggers_low_res_full_sky',
                    'HERA331_and_Joshs_outriggers_and_3_inriggers_low_res_full_sky']


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

                    
#%% Load results of calcuations
scriptDirectory = os.path.dirname(os.path.abspath(__file__))
allSpecs = []
allPSFs = []
allmapNoiseCovariances = []
allAtransNinvAs = []
allnoiseEvals = []
allnoiseEvecs = []
allCoords = []
allOmniAtransA = []
allGainVariances = []
allNZeroEVs = []

for result in ResultsToExamine:
    print "Now loading from " + scriptDirectory + '/Results/' + result + '/'    
    resultsDirectory = scriptDirectory + '/Results/' + result + '/'
    s, times, ps, Dmatrix, PSF, coaddedMap, pointSourcePSF = MapMats.loadAllResults(resultsDirectory)
    allSpecs.append(s)
    allPSFs.append(PSF)
    allmapNoiseCovariances.append(np.load(resultsDirectory + "mapNoiseCovariance.npy"))
    allAtransNinvAs.append(np.load(resultsDirectory + "AtransNinvA.npy"))
    allnoiseEvals.append(np.load(resultsDirectory + "noiseEvals.npy"))
    allnoiseEvecs.append(np.load(resultsDirectory + "noiseEvecs.npy"))
    allCoords.append(pickle.load(open(s.resultsFolder + "coords.p","rb")))
    allOmniAtransA.append(np.load(resultsDirectory + "omniAtransA.npy"))
    allGainVariances.append(np.load(resultsDirectory + "gainVariances.npy"))
    allNZeroEVs.append(str(np.loadtxt(resultsDirectory + "nZeroEVs.txt")))

    

#%% Plotting

plt.close('all')
for i,folder in enumerate(ResultsToExamine):
    s = allSpecs[i]
    
    plt.figure(figsize=(14,6))
    plt.subplot(1, 2, 1)    

    #plt.imshow(AtransA, interpolation='none')
    #plt.imshow(np.linalg.pinv(AtransA)[0:s.nAntennas,0:s.nAntennas], interpolation='none')   
  #  plt.scatter(s.antennaPositions[:,0], s.antennaPositions[:,1], c=gainVariances, s=100, norm=matplotlib.colors.LogNorm())
    plt.scatter(s.antennaPositions[:,0], s.antennaPositions[:,1], c=allGainVariances[i], s=100)
    plt.colorbar()
    print "AtransA has " + allNZeroEVs[i] + " zero eigenvalues."
    datacursor(display='single',formatter="x={x:.2f}\ny={y:.2f}".format)

    plt.subplot(1, 2, 2)
    plt.semilogy(np.real(allnoiseEvals[i]))
    plt.xlabel('Eigenvalue Number',fontsize=16)
    plt.ylabel('Eigenvalue of $\mathbf{A}^\dagger\mathbf{N}^{-1}\mathbf{A}$',fontsize=16)
    #datacursor(display='single',formatter="N={x:.2f}\nEV={y:.2e}".format)

    
