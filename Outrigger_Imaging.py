import numpy as np
import healpy as hp
import math
import sys,os,os.path
sys.path.append('/Users/jsdillon/Desktop/')
import matplotlib
import matplotlib.pyplot as plt
from Joint_Mapmaking_Power_Spectrum_Pipeline.Source.Specifications import Specifications
from Joint_Mapmaking_Power_Spectrum_Pipeline.Source import Geometry
from Joint_Mapmaking_Power_Spectrum_Pipeline.Mapmaker import Mapmaker
import scipy.constants as const

#sys.stdout = sys.__stdout__    
print sys.stdout
plt.close('all')

print "testing"


scriptDirectory = os.path.dirname(os.path.abspath(__file__))

Mapmaker(configFile = "outrigger_config.txt", mainDirectory = scriptDirectory)

