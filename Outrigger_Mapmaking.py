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
import smtplib
from email.mime.text import MIMEText
import email.Utils

scriptDirectory = os.path.dirname(os.path.abspath(__file__))
os.chdir(scriptDirectory)

def SortedEigensystem(matrix):
    """Returns the eigensystem of the input matrix where eigenvalues and eigenvectors are sorted by descending absolute value."""
    evals,evecs = np.linalg.eig(matrix)
    indices = np.argsort(np.abs(evals))[::-1]   
    return evals[indices], evecs[:,indices]

def saveQuantitiesForArrayComparison(resultsDirectory):
    """Computes and saves interesting quantities for the outrigger comparsion in the same folder as the mapmaking results."""
    print "Now computing other interesting mapmaking and array quantities..."
    s, times, ps, Dmatrix, PSF, coaddedMap, pointSourcePSF = MapMats.loadAllResults(resultsDirectory)

    #Info about calibratability    
    redundantInfo = oi.RedundantInfo()
    reds = oai.compute_reds(s.antennaPositions)
    redundantInfo.init_from_reds(reds, s.antennaPositions)
    AtransA = redundantInfo.At.dot(redundantInfo.At.T).toarray()
    BtransB = redundantInfo.Bt.dot(redundantInfo.Bt.T).toarray()    
    np.save(resultsDirectory + "omniAtransA",AtransA)
    np.save(resultsDirectory + "omniBtransB",BtransB)
    gainVariances = np.diag(np.linalg.pinv(AtransA)[0:s.nAntennas,0:s.nAntennas])
    np.save(resultsDirectory + "gainVariances",gainVariances)
    np.savetxt(resultsDirectory + "nZeroEVs.txt", [np.sum(np.abs(np.linalg.eigvals(XtransX)<1e-10)) for XtransX in [AtransA, BtransB]], fmt='%.d')    
    
    coords = Geometry.Coordinates(s)
    if s.useAdaptiveHEALPixForPSF: coords.convertToAdaptiveHEALPix(s, times)
    print "There are " + str(coords.nPSFPixels) + " pixels in the PSF and " + str(coords.nFacetPixels) + " in the map."
    pickle.dump(coords, open(resultsDirectory + "coords.p","wb"), protocol = -1)    
    mapNoiseCovariance = np.dot(PSF[:,coords.facetIndexLocationsInPSFIndexList],np.transpose(np.diag(Dmatrix)))
    AtransNinvA = np.dot(np.diag(Dmatrix**-1),PSF[:,coords.facetIndexLocationsInPSFIndexList])
    posVarIndicies = (np.diag(AtransNinvA)>0)
    AtransNinvA = AtransNinvA[posVarIndicies,:][:,posVarIndicies] #Deal with the fact that if I map the whole sky, some pixels are not sampled.
    mapNoiseCovariance = mapNoiseCovariance[posVarIndicies,:][:,posVarIndicies]
    evals, evecs = SortedEigensystem(AtransNinvA)
    np.save(resultsDirectory + "mapNoiseCovariance",mapNoiseCovariance)
    np.save(resultsDirectory + "AtransNinvA", AtransNinvA)
    np.save(resultsDirectory + "noiseEvals", evals)
    np.save(resultsDirectory + "noiseEvecs", evecs)
    np.save(resultsDirectory + "posVarIndicies", posVarIndicies)    
    print "Finished with " + resultsDirectory + "\n\n\n"



highResMapNSIDE = 256
allTests = [
#%% Full Sky Single Integration
            {'title': 'HERA-127 Full Sky', 'folder': 'HERA127_FullSky/', 'hexNum': 7,
                 'mapNSIDE': 64, 'useAdaptiveHEALPixForPSF': True, 'facetSize': 10},
                 
            {'title': 'HERA-331 Full Sky', 'folder': 'HERA331_FullSky/', 
                 'mapNSIDE': 64, 'useAdaptiveHEALPixForPSF': True, 'facetSize': 10},
                 
            {'title': 'HERA-331 + Fiducial Inriggers Full Sky', 'folder': 'HERA331_FidInriggers_FullSky/', 
                 'fiducialInriggers': True, 'mapNSIDE': 64, 'useAdaptiveHEALPixForPSF': True, 'facetSize': 10},
                 
#            {'title': 'HERA-331 + 3 Corner Inrigger Pairs Full Sky', 'folder': 'HERA331_3CornerInriggers_FullSky/', 
#                 'halfCornerInriggers': True, 'mapNSIDE': 64, 'useAdaptiveHEALPixForPSF': True, 'facetSize': 10},
                 
            {'title': 'HERA-331 + 6 Corner Inrigger Pairs Full Sky', 'folder': 'HERA331_6CornerInriggers_FullSky/', 
                 'fullCornerInriggers': True, 'mapNSIDE': 64, 'useAdaptiveHEALPixForPSF': True, 'facetSize': 10},
                 
            {'title': 'HERA-331 Split Core Full Sky', 'folder': 'HERA331_SplitCore_FullSky/', 'SplitCore': True, 'mapNSIDE': 64, 
                 'useAdaptiveHEALPixForPSF': True, 'facetSize': 10},
                 
            {'title': 'HERA-331 Nine-Way Split Core Full Sky', 'folder': 'HERA331_NineWaySplitCore_FullSky/', 'NineWaySplitCore': True, 'mapNSIDE': 64, 
                 'useAdaptiveHEALPixForPSF': True, 'facetSize': 10},
                 
#%% Full Sky Multiple Integrations
            {'title': 'HERA-331 Full Sky Long Integration', 'folder': 'HERA331_FullSky_LongInt/', 
                 'MaximumAllowedAngleFromFacetCenterToPointingCenter': 5, 'mapNSIDE': 64, 'useAdaptiveHEALPixForPSF': True, 'facetSize': 10},
                 
            {'title': 'HERA-331 + Fiducial Inriggers Full Sky Long Integration', 'folder': 'HERA331_FidInriggers_FullSky_LongInt/', 
                 'MaximumAllowedAngleFromFacetCenterToPointingCenter': 5, 'fiducialInriggers': True, 'mapNSIDE': 64, 'useAdaptiveHEALPixForPSF': True, 'facetSize': 10},
                 
#            {'title': 'HERA-331 + 3 Corner Inrigger Pairs Full Sky Long Integration', 'folder': 'HERA331_3CornerInriggers_FullSky_LongInt/', 
#                 'MaximumAllowedAngleFromFacetCenterToPointingCenter': 5, 'halfCornerInriggers': True, 'mapNSIDE': 64, 'useAdaptiveHEALPixForPSF': True, 'facetSize': 10},
                 
            {'title': 'HERA-331 + 6 Corner Inrigger Pairs Full Sky Long Integration', 'folder': 'HERA331_6CornerInriggers_FullSky_LongInt/', 
                 'MaximumAllowedAngleFromFacetCenterToPointingCenter': 5, 'fullCornerInriggers': True, 'mapNSIDE': 64, 'useAdaptiveHEALPixForPSF': True, 'facetSize': 10},
                 
            {'title': 'HERA-331 Split Core Full Sky Long Integration', 'folder': 'HERA331_SplitCore_FullSky_LongInt/', 
                 'MaximumAllowedAngleFromFacetCenterToPointingCenter': 5, 'SplitCore': True, 'mapNSIDE': 64, 'useAdaptiveHEALPixForPSF': True, 'facetSize': 10},
                 
            {'title': 'HERA-331 Nine-Way Split Core Full Sky Long Integration', 'folder': 'HERA331_NineWaySplitCore_FullSky_LongInt/', 
                 'MaximumAllowedAngleFromFacetCenterToPointingCenter': 5, 'NineWaySplitCore': True, 'mapNSIDE': 64, 'useAdaptiveHEALPixForPSF': True, 'facetSize': 10},

#%% High Res Partial Sky
            {'title': 'HERA-331 High-Res', 'folder': 'HERA331_HighRes/', 
                 'mapNSIDE': highResMapNSIDE, 'useAdaptiveHEALPixForPSF': False, 'facetSize': 10, 'PSFextensionBeyondFacetFactor': 2},
                         
#            {'title': 'HERA-331 + Fidicuial Outriggers High-Res', 'folder': 'HERA331_FidOutriggers_HighRes/', 
#                 'LoadHERAOutriggers': True, 'mapNSIDE': highResMapNSIDE, 'useAdaptiveHEALPixForPSF': False, 'facetSize': 10, 'PSFextensionBeyondFacetFactor': 2},
                 
#            {'title': 'HERA-331 + Redundant Outriggers High-Res', 'folder': 'HERA331_RedOutriggers_HighRes/', 
#                 'RedundantOutriggers': True, 'mapNSIDE': highResMapNSIDE, 'useAdaptiveHEALPixForPSF': False, 'facetSize': 10, 'PSFextensionBeyondFacetFactor': 2},
                 
            {'title': 'HERA-331 + Fidicuial Outriggers + Fidicuial Inriggers High-Res', 'folder': 'HERA331_FidOutriggers_FidInriggers_HighRes/', 
                 'LoadHERAOutriggers': True, 'fiducialInriggers': True, 'mapNSIDE': highResMapNSIDE, 'useAdaptiveHEALPixForPSF': False, 'facetSize': 10, 'PSFextensionBeyondFacetFactor': 2},
                 
#            {'title': 'HERA-331 + Redundant Outriggers + 3 Corner Inrigger Pairs High-Res', 'folder': 'HERA331_RedOutriggers_3PairInriggers_HighRes/', 
#                 'RedundantOutriggers': True, 'halfCornerInriggers': True, 'mapNSIDE': highResMapNSIDE, 'useAdaptiveHEALPixForPSF': False, 'facetSize': 10, 'PSFextensionBeyondFacetFactor': 2},
                 
            {'title': 'HERA-331 + Redundant Outriggers + 6 Corner Inrigger Pairs High-Res', 'folder': 'HERA331_RedOutriggers_6PairInriggers_HighRes/', 
                 'RedundantOutriggers': True, 'fullCornerInriggers': True, 'mapNSIDE': highResMapNSIDE, 'useAdaptiveHEALPixForPSF': False, 'facetSize': 10, 'PSFextensionBeyondFacetFactor': 2},
                 
            {'title': 'HERA-331 Split Core + Split Outriggers High-Res', 'folder': 'HERA331_SplitCore_SplitOutriggers_HighRes/', 
                 'SplitCore': True, 'SplitCoreOutriggers': True, 'mapNSIDE': highResMapNSIDE, 'useAdaptiveHEALPixForPSF': False, 'facetSize': 10, 'PSFextensionBeyondFacetFactor': 2},

            {'title': 'HERA-331 Nine-Way Split Core + Split Outriggers High-Res', 'folder': 'HERA331_NineWaySplitCore_SplitOutriggers_HighRes/', 
                 'NineWaySplitCore': True, 'SplitCoreOutriggers': True, 'mapNSIDE': highResMapNSIDE, 'useAdaptiveHEALPixForPSF': False, 'facetSize': 10, 'PSFextensionBeyondFacetFactor': 2},

#%% High Res Partial Sky Long Int
            {'title': 'HERA-331 High-Res Long Integration', 'folder': 'HERA331_HighRes_LongInt/', 
                 'MaximumAllowedAngleFromFacetCenterToPointingCenter': 5, 'mapNSIDE': highResMapNSIDE, 'useAdaptiveHEALPixForPSF': False, 'facetSize': 10, 'PSFextensionBeyondFacetFactor': 2},
                         
#            {'title': 'HERA-331 + Fidicuial Outriggers High-Res Long Integration', 'folder': 'HERA331_FidOutriggers_HighRes_LongInt/', 
#                 'MaximumAllowedAngleFromFacetCenterToPointingCenter': 5, 'LoadHERAOutriggers': True, 'mapNSIDE': highResMapNSIDE, 'useAdaptiveHEALPixForPSF': False, 'facetSize': 10, 'PSFextensionBeyondFacetFactor': 2},
                 
#            {'title': 'HERA-331 + Redundant Outriggers High-Res Long Integration', 'folder': 'HERA331_RedOutriggers_HighRes_LongInt/', 
#                 'MaximumAllowedAngleFromFacetCenterToPointingCenter': 5, 'RedundantOutriggers': True, 'mapNSIDE': highResMapNSIDE, 'useAdaptiveHEALPixForPSF': False, 'facetSize': 10, 'PSFextensionBeyondFacetFactor': 2},
                 
            {'title': 'HERA-331 + Fidicuial Outriggers + Fidicuial Inriggers High-Res Long Integration', 'folder': 'HERA331_FidOutriggers_FidInriggers_HighRes_LongInt/', 
                 'MaximumAllowedAngleFromFacetCenterToPointingCenter': 5, 'LoadHERAOutriggers': True, 'fiducialInriggers': True, 'mapNSIDE': highResMapNSIDE, 'useAdaptiveHEALPixForPSF': False, 'facetSize': 10, 'PSFextensionBeyondFacetFactor': 2},
                 
#            {'title': 'HERA-331 + Redundant Outriggers + 3 Corner Inrigger Pairs High-Res Long Integration', 'folder': 'HERA331_RedOutriggers_3PairInriggers_HighRes_LongInt/', 
#                 'MaximumAllowedAngleFromFacetCenterToPointingCenter': 5, 'RedundantOutriggers': True, 'halfCornerInriggers': True, 'mapNSIDE': highResMapNSIDE, 'useAdaptiveHEALPixForPSF': False, 'facetSize': 10, 'PSFextensionBeyondFacetFactor': 2},
                 
            {'title': 'HERA-331 + Redundant Outriggers + 6 Corner Inrigger Pairs High-Res Long Integration', 'folder': 'HERA331_RedOutriggers_6PairInriggers_HighRes_LongInt/', 
                 'MaximumAllowedAngleFromFacetCenterToPointingCenter': 5, 'RedundantOutriggers': True, 'fullCornerInriggers': True, 'mapNSIDE': highResMapNSIDE, 'useAdaptiveHEALPixForPSF': False, 'facetSize': 10, 'PSFextensionBeyondFacetFactor': 2},
                 
            {'title': 'HERA-331 Split Core + Split Outriggers High-Res Long Integration', 'folder': 'HERA331_SplitCore_SplitOutriggers_HighRes_LongInt/', 
                 'MaximumAllowedAngleFromFacetCenterToPointingCenter': 5, 'SplitCore': True, 'SplitCoreOutriggers': True, 'mapNSIDE': highResMapNSIDE, 'useAdaptiveHEALPixForPSF': False, 'facetSize': 10, 'PSFextensionBeyondFacetFactor': 2},

            {'title': 'HERA-331 Nine-Way Split Core + Split Outriggers High-Res Long Integration', 'folder': 'HERA331_NineWaySplitCore_SplitOutriggers_HighRes_LongInt/', 
                 'MaximumAllowedAngleFromFacetCenterToPointingCenter': 5, 'NineWaySplitCore': True, 'SplitCoreOutriggers': True, 'mapNSIDE': highResMapNSIDE, 'useAdaptiveHEALPixForPSF': False, 'facetSize': 10, 'PSFextensionBeyondFacetFactor': 2}]

            

#TODO: look at half-wavelength dithering specifically in the N-S direction...or extensions to the "paired" inriggers
pickle.dump(allTests, open(scriptDirectory + "/Results/allTests.p","wb"))


def emailMe(test):
    fromaddr = 'joshspythonnotifications@gmail.com'
    toaddrs  = 'jsdillon+python@gmail.com'
    message = str(test)
    msg = MIMEText(message)
    msg['Subject'] = 'Finished test: ' + test['title']
    msg['Message-id'] = email.Utils.make_msgid()
    msg['From'] = fromaddr
    msg['To'] = toaddrs
    username = 'joshspythonnotifications'
    password = 'TqU2rLmjPpme9uTHLrn'
    server = smtplib.SMTP('smtp.gmail.com:587')
    server.starttls()
    server.login(username,password)
    server.sendmail(fromaddr, toaddrs, msg.as_string())
    server.quit()


def Outrigger_Mapmaking(testCase = None):
    """This function will run either one array configuration, or all of them if testCase == None."""
    
    if len(sys.argv) > 1: testCase = int(sys.argv[1])
    if testCase is None: print "Now working on all array configurations.\n\n"
    
    for case, test in enumerate(allTests):
        if case == testCase or testCase is None:
            HexArray(**test)
            SampleObservationDataGenerator(configFile = "outrigger_config.txt")
            resultsDirectory = Mapmaker(resultsFolder = scriptDirectory + "/Results/" + test['folder'], configFile = "outrigger_config.txt", mainDirectory = scriptDirectory, **test)
            saveQuantitiesForArrayComparison(resultsDirectory)
            emailMe(test)


#%%
if __name__ == "__main__":
    Outrigger_Mapmaking(testCase = 17)    
    Outrigger_Mapmaking(testCase = 18)
    Outrigger_Mapmaking(testCase = 19)
    Outrigger_Mapmaking(testCase = 20)
#    Outrigger_Mapmaking(testCase = 6)
#    Outrigger_Mapmaking(testCase = 18)
#    Outrigger_Mapmaking(testCase = 19)
  #  Outrigger_Mapmaking()
            
#    plt.hist(
    
    #Future Plans for this:
    #Examine the effect of inriggers on low-resolution
    #Examine the effect of outriggers on high resolution
        #Look at different dithering schemes
    #Examine metrics for calibratability
    #Examine array rotation
    
    #calibratibility
    
    #invert polarization beams
    
    #Adrian's Idea for the Paper:
    #Intuition/metrics/rules of thumb
        #grating lobes, 
        #redundant calibratibility
        #high resolution
        #beam symmetry for uniform weighted beam
    

# How sample variance limited are we and how much do we buy with the extra modes
# How good is good enough on the grating lobes
# 