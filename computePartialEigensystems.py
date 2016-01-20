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

scriptDirectory = os.path.dirname(os.path.abspath(__file__))
allTests = pickle.load(open(scriptDirectory + "/Results/allTests.p","rb"))

rangesOfInterest = [np.arange(0,100), 
                    np.arange(0,630), 
                    np.arange(0,1000),
                    np.arange(0,1200),
                    np.arange(1200,5900),
                    np.arange(1300,1500),
                    np.arange(630,3000),
                    np.arange(1000,2000), 
                    np.arange(1000,3000),
                    np.arange(1000,4000),
                    np.arange(3000,4000),
                    np.arange(4000,5900),
                    np.arange(4150,4750),
                    np.arange(4750,5900)]

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

#%% Compute and save eigenvalue segments
        result = allResults[test['title']]
        for thisRange in rangesOfInterest:
            resultsFile =  scriptDirectory + '/Results/' + test['folder'] + '/partialEigenspace_' + str(thisRange[0]) + "-" + str(thisRange[-1])
            if not os.path.isfile(resultsFile):
                print "Working on " + test['title'] + ' for eigenvalues ' + str(thisRange[0]) + "-" + str(thisRange[-1])
                partialEigenspace = np.einsum('ij,jj,jk->ik',result['noiseEvecs'][:,thisRange],np.diag(result['noiseEvals'][thisRange]),result['noiseEvecs'][:,thisRange].conj().T) 
                np.save(resultsFile, np.diag(partialEigenspace))

        del allResults[test['title']]
