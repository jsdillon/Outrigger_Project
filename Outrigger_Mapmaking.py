import numpy as np
import healpy as hp
import math
import sys,os,os.path
import matplotlib
import matplotlib.pyplot as plt
from MAPS21cm.Specifications import Specifications
from MAPS21cm import Geometry
from MAPS21cm.Mapmaker import Mapmaker
import MAPS21cm.MatricesForMapmaking as MapMats
import scipy.constants as const
from InstrumentData.HexArray import HexArray
from ObservationData.SampleObservationDataGenerator import SampleObservationDataGenerator
import omnical.info as oi
import omnical.arrayinfo as oai
import scipy.sparse
import cPickle as pickle

scriptDirectory = os.path.dirname(os.path.abspath(__file__))

def SortedEigensystem(matrix):
    """Returns the eigensystem of the input matrix where eigenvalues and eigenvectors are sorted by descending absolute value."""
    evals,evecs = np.linalg.eig(matrix)
    indices = np.argsort(np.abs(evals))[::-1]   
    return evals[indices], evecs[:,indices]

def saveQuantitiesForArrayComparison(resultsDirectory):
    """Computes and saves interesting quantities for the outrigger comparsion in the same folder as the mapmaking results."""
    print "Now computing other interesting mapmaking and array quantities..."
    s, times, ps, Dmatrix, PSF, coaddedMap, pointSourcePSF = MapMats.loadAllResults(resultsDirectory)

    redundantInfo = oi.RedundantInfo()
    reds = oai.compute_reds(s.antennaPositions)
    redundantInfo.init_from_reds(reds, s.antennaPositions)
    AtransA = redundantInfo.At.dot(redundantInfo.At.T).toarray()
    np.save(resultsDirectory + "omniAtransA",AtransA)
    gainVariances = np.diag(np.linalg.pinv(AtransA)[0:s.nAntennas,0:s.nAntennas])
    gainVariances = gainVariances/gainVariances.min()
    np.save(resultsDirectory + "gainVariances",gainVariances)
    np.savetxt(resultsDirectory + "nZeroEVs.txt", [len(AtransA) - np.linalg.matrix_rank(AtransA, tol=1e-10)], fmt='%.d')
    
    coords = Geometry.Coordinates(s)
    pickle.dump(coords, open(resultsDirectory + "coords.p","wb"))    
    mapNoiseCovariance = np.dot(PSF[:,coords.mapIndexLocationsInExtendedIndexList],np.transpose(np.diag(Dmatrix)))
    AtransNinvA = np.dot(np.diag(Dmatrix**-1),PSF)    
    evals, evecs = SortedEigensystem(AtransNinvA)
    np.save(resultsDirectory + "mapNoiseCovariance",mapNoiseCovariance)
    np.save(resultsDirectory + "AtransNinvA", AtransNinvA)
    np.save(resultsDirectory + "noiseEvals", evals)
    np.save(resultsDirectory + "noiseEvecs", evecs)
    print "Finished with " + resultsDirectory + "\n\n\n"
    
def Outrigger_Mapmaking(testCase = None):
    """This function will run either one array configuration, or all of them if testCase == None."""
    
    #%% Basic HERA 331
    if testCase == 0 or testCase == None:
        HexArray()
        SampleObservationDataGenerator(configFile = "outrigger_config.txt")
        resultsDirectory = Mapmaker(resultsFolder = scriptDirectory + "/Results/HERA331/", configFile = "outrigger_config.txt", mainDirectory = scriptDirectory)
        saveQuantitiesForArrayComparison(resultsDirectory)
    
    
    #%% HERA 331 + 3 inriggers
    if testCase == 1 or testCase == None:
        HexArray(ditheredInriggers=True)
        SampleObservationDataGenerator(configFile = "outrigger_config.txt")
        resultsDirectory = Mapmaker(resultsFolder = scriptDirectory + "/Results/HERA331_and_3_inriggers/", configFile = "outrigger_config.txt", mainDirectory = scriptDirectory)
        saveQuantitiesForArrayComparison(resultsDirectory)
    
    
    #%% HERA 331 with my outrigger design
    if testCase == 2 or testCase == None:
        HexArray(JoshsOutriggers=True)
        SampleObservationDataGenerator(configFile = "outrigger_config.txt")
        resultsDirectory = Mapmaker(resultsFolder = scriptDirectory + "/Results/HERA331_and_Joshs_outriggers/", configFile = "outrigger_config.txt", mainDirectory = scriptDirectory)
        saveQuantitiesForArrayComparison(resultsDirectory)
        
    
    #%% HERA 331 with my outrigger design and 3 inriggers
    if testCase == 3 or testCase == None:
        HexArray(JoshsOutriggers=True,ditheredInriggers=True)
        SampleObservationDataGenerator(configFile = "outrigger_config.txt")
        resultsDirectory = Mapmaker(resultsFolder = scriptDirectory + "/Results/HERA331_and_Joshs_outriggers_and_3_inriggers/", configFile = "outrigger_config.txt", mainDirectory = scriptDirectory)
        saveQuantitiesForArrayComparison(resultsDirectory)
    
    
    #%% HERA 331 with full sky at NSIDE=16
    if testCase == 4 or testCase == None:    
        HexArray()
        SampleObservationDataGenerator(configFile = "outrigger_config.txt")
        resultsDirectory = Mapmaker(resultsFolder = scriptDirectory + "/Results/HERA331_low_res_full_sky/", configFile = "outrigger_config.txt", mainDirectory = scriptDirectory, 
                                    facetSize = 180, mapNSIDE = 16, PSFextensionBeyondFacetFactor = 1)
        saveQuantitiesForArrayComparison(resultsDirectory)
        
    
    #%% HERA 331 with 3 inriggers and full sky at NSIDE=16
    if testCase == 5 or testCase == None:    
        HexArray(ditheredInriggers=True)
        SampleObservationDataGenerator(configFile = "outrigger_config.txt")
        resultsDirectory = Mapmaker(resultsFolder = scriptDirectory + "/Results/HERA331_and_3_inriggers_low_res_full_sky/", configFile = "outrigger_config.txt", mainDirectory = scriptDirectory, 
                                    facetSize = 180, mapNSIDE = 16, PSFextensionBeyondFacetFactor = 1)                            
        saveQuantitiesForArrayComparison(resultsDirectory)                                                        
        
    
    #%% HERA 331 with 3 inriggers, my outrigger design, and full sky at NSIDE=16
    if testCase == 6 or testCase == None:    
        HexArray(JoshsOutriggers=True,ditheredInriggers=True)
        SampleObservationDataGenerator(configFile = "outrigger_config.txt")
        resultsDirectory = Mapmaker(resultsFolder = scriptDirectory + "/Results/HERA331_and_Joshs_outriggers_and_3_inriggers_low_res_full_sky/", configFile = "outrigger_config.txt", mainDirectory = scriptDirectory, 
                                    facetSize = 180, mapNSIDE = 16, PSFextensionBeyondFacetFactor = 1)                            
        saveQuantitiesForArrayComparison(resultsDirectory)          


#%%
if __name__ == "__main__":
    Outrigger_Mapmaking()

    
    #Future Plans for this:
    #Examine the effect of inriggers on low-resolution
    #Examine the effect of outriggers on high resolution
        #Look at different dithering schemes
    #Examine metrics for calibratability
    #Examine array rotation
    
    #Adrian's Idea for the Paper:
    #Intuition/metrics/rules of thumb
        #grating lobes, 
        #redundant calibratibility
        #high resolution
        #beam symmetry for uniform weighted beam
    
