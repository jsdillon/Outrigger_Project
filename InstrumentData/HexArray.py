# TEST CODE FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon
#
# This code generates a file of antenna positions and a files with baselines, baseline info, and redundancies

import os
import numpy as np
import cPickle as pickle
import matplotlib
import matplotlib.pyplot as plt
#from mpldatacursor import datacursor
plt.close('all')

def HexArray(Separation = 14.6, hexNum = 11, 
             ditheredInriggers = False, 
             AaronsInriggers = False, 
             JoshsOutriggers = False, 
             AaronsOutriggers = False, 
             LoadHERAOutriggers = False, 
             redundantTriangleInriggers = False,
             redundantPairInriggers = False,
             redundantHexInriggers = False):
    
    precisionFactor = 1000000        
    
    #Calculating Positions:
    #Main Hex
    positions = [];
    for row in range(hexNum-1,-(hexNum),-1):
        for col in range(2*hexNum-abs(row)-1):
            xPos = ((-(2*hexNum-abs(row))+2)/2.0 + col)*Separation;
            yPos = row*Separation*3**.5/2;
            positions.append([xPos, yPos, 0])
    nCore = len(positions)

    right = Separation*np.asarray([1,0,0])
    up = Separation*np.asarray([0,1,0])
    upRight = Separation*np.asarray([.5,3**.5/2,0])
    upLeft = Separation*np.asarray([-.5,3**.5/2,0])

    #Inriggers
    if ditheredInriggers:
        print "adding Josh's inriggers"        
        if hexNum % 2 == 1:
            positions.append(.5*right + hexNum * (upLeft + upRight)/2.0 + .5 * upLeft)
            positions.append(hexNum*right + (hexNum-1)/2*upLeft + .5*upRight)
            positions.append(-hexNum*right + (hexNum+1)/2*upRight  - .5*right)
        else:
            positions.append(hexNum * (upLeft + upRight)/2.0 + .5 * upLeft)
            positions.append(hexNum*right + (hexNum)/2*upLeft + .5*upRight)
            positions.append(-hexNum*right + (hexNum)/2*upRight  - .5*right)

    #Aaron's alternate inriggers
    if AaronsInriggers:
        print "adding Aaron's inriggers"        
        if hexNum % 2 == 1:
            positions.append(upRight + ((hexNum-1)/2 + 1.0/3**.5/2.0) * (upLeft + upRight))
            positions.append(right + ((hexNum-1)/2 + 1.0/3**.5/2.0) * (right + upRight))
            positions.append(right + ((hexNum+1)/2 + 1.0/3**.5/2.0) * (-right + upLeft))
        else:
            positions.append(((hexNum)/2 + 1.0/3**.5/2.0) * (upLeft + upRight))
            positions.append(((hexNum)/2 + 1.0/3**.5/2.0) * (right + upRight))
            positions.append(((hexNum)/2 + 1.0/3**.5/2.0) * (-right + upLeft))

    #Redundant Inriggers
    if redundantTriangleInriggers or redundantPairInriggers or redundantHexInriggers:
        hexNumInrigger = 2
        if hexNum % 2 == 1:            
            inriggerCenters = [((hexNum-1)/2) * (upLeft + upRight) + upLeft + upLeft + (right + upRight)/3,       #top
                               -((hexNum-1)/2) * (upLeft + upRight) - upLeft - upLeft - 2*(right + upRight)/3,    #bottom
                                -(-(hexNum-1)/2 * (upRight + right) - upRight - right - (right + upRight)/3),     #top right  
                                -((hexNum-1)/2 * (upRight + right) + upRight  + 2*(right + upRight)/3),      #bottom left
                                (hexNum-1)/2 * (upLeft - right) + upLeft + upRight - right - 2*(right + upRight)/3,   #top left
                                -(hexNum-1)/2 * (upLeft - right) - upLeft + right + (right + upRight)/3]        #top right
        else:
            inriggerCenters = [hexNum/2 * (upLeft + upRight) + upLeft  + (right + upRight)/3,       #top
                               -hexNum/2 * (upLeft + upRight) - upLeft - 2*(right + upRight)/3,    #bottom
                                hexNum/2 * (upRight + right)   + 2*(right + upRight)/3,      #top right
                                -hexNum/2 * (upRight + right) - right - (right + upRight)/3,       #bottom left
                                hexNum/2 * (upLeft - right) + upRight - right - 2*(right + upRight)/3,   #top left
                                -hexNum/2 * (upLeft - right) - upLeft  + (right + upRight)/3]        #top right
        for newCenter in inriggerCenters:            
            for row in range(hexNumInrigger-1,-(hexNumInrigger),-1):
                for col in range(2*hexNumInrigger-abs(row)-1):
                    xPos = ((-(2*hexNumInrigger-abs(row))+2)/2.0 + col)*Separation + newCenter[0]
                    yPos = row*Separation*3**.5/2 + newCenter[1]
                    if redundantHexInriggers: positions.append([xPos, yPos, 0])
                    elif redundantTriangleInriggers:
                        if (xPos**2+yPos**2)**.5 < np.linalg.norm(newCenter) or np.array_equal(np.array([xPos, yPos, 0]), newCenter): positions.append([xPos, yPos, 0])
                    elif redundantPairInriggers:
                        if (xPos**2+yPos**2)**.5 < np.linalg.norm(newCenter): positions.append([xPos, yPos, 0])
                                                    

    
    #Outriggers
    if JoshsOutriggers:
        outriggerHexNum = 3
        for row in range(outriggerHexNum-1,-(outriggerHexNum),-1):
            for col in range(2*outriggerHexNum-abs(row)-1):
                xPos = ((-(2*outriggerHexNum-abs(row))+2)/2.0 + col)*Separation*(np.floor(1.5*hexNum));
                yPos = row*Separation*(np.floor(1.5*hexNum))*3**.5/2;
                if xPos != 0 or yPos != 0:    
                    positions.append([xPos, yPos, 0])

    if AaronsOutriggers:
        outriggerHexNum = 3
        for row in range(outriggerHexNum-1,-(outriggerHexNum),-1):
            for col in range(2*outriggerHexNum-abs(row)-1):
                yPos = ((-(2*outriggerHexNum-abs(row))+2)/2.0 + col)*Separation*(3**.5*(hexNum-1));
                xPos = row*Separation*(3**.5*(hexNum-1))*3**.5/2;
                if xPos != 0 or yPos != 0:    
                    positions.append([xPos, yPos, 0])
        
    if LoadHERAOutriggers:
        outriggerPositions = np.loadtxt('hexconfig352_outriggers_only.dat')
        for outrigger in outriggerPositions:
            positions.append(np.append(outrigger,0))
        precisionFactor = 10000

    #Calculating Baselines    
   
    nAntennas = len(positions)
    baselines = []
    baselinePairs = []
    #print "WARNING: DOUBLING UNIQUE BASELINES!!!"
    for ant1 in range(nAntennas):
        startAnt2 = ant1+1
        if __name__ == "__main__":
            print "WARNING: DOUBLING UNIQUE BASELINES!!!"
            startAnt2 = 0
        for ant2 in range(startAnt2,nAntennas):
            baselines.append((int(np.round(precisionFactor*(positions[ant1][0]-positions[ant2][0]))), int(np.round(precisionFactor*(positions[ant1][1]-positions[ant2][1]))), int(np.round(precisionFactor*(positions[ant1][2]-positions[ant2][2])))))
            baselinePairs.append((ant1, ant2))
    
    #Calculating Unique Baselines
    baselineDict = {}
    for b in range(len(baselines)):
        if baselineDict.has_key(baselines[b]):
            baselineDict[baselines[b]].append(baselinePairs[b])
        else:
            baselineDict[baselines[b]] = [baselinePairs[b]]
    
    print "With", len(positions), "antennas there are", len(baselineDict.items()), "unique baselines."
    
    #Saving results
    scriptDirectory = os.path.dirname(os.path.abspath(__file__))
    np.savetxt(scriptDirectory + "/antenna_positions.dat",np.asarray(positions))
    np.savetxt(scriptDirectory + "/all_baselines.dat",np.asarray(baselines)/(1.0*precisionFactor))
    np.savetxt(scriptDirectory + "/all_baseline_pairs.dat",np.asarray(baselinePairs),fmt='%i')
    np.savetxt(scriptDirectory + "/unique_baselines.dat",np.asarray([uniqueBaseline[0] for uniqueBaseline in baselineDict.items()])/(1.0*precisionFactor))
    np.savetxt(scriptDirectory + "/redundancy.dat",np.asarray([len(uniqueBaseline[1]) for uniqueBaseline in baselineDict.items()]),fmt='%i')
    
    antennaPairDict = {}
    for uniqueIndex in range(len(baselineDict.items())):
        allItems = baselineDict.items()
        for antennaPair in allItems[uniqueIndex][1]:
            antennaPairDict[antennaPair] = uniqueIndex
    pickle.dump(antennaPairDict, open(scriptDirectory + "/antennaPairUniqueBaselineIndexDict.p", 'w'))
    
    if __name__ == "__main__":
#        fig = plt.gcf()
#        plt.ylim(-600,600)
#        plt.xlim(-600,600)
#        for antPos in positions:
#            fig.gca().add_artist(plt.Circle((antPos[0],antPos[1]),7,fc='0.3',ec='k'))
        #datacursor(display='single')
        
        uniqueBaselines = np.asarray([uniqueBaseline[0] for uniqueBaseline in baselineDict.items()])/(1.0*precisionFactor)
        redundancy = np.asarray([len(uniqueBaseline[1]) for uniqueBaseline in baselineDict.items()])
#        fig2, axes = plt.subplots(1,2)
        plt.figure()
        plt.scatter(np.asarray(positions)[:,0]/Separation,np.asarray(positions)[:,1]/Separation)        
        #datacursor(display='single')
        plt.figure()        
        plt.scatter(uniqueBaselines[:,0]/1.0/Separation, uniqueBaselines[:,1]/1.0/Separation,c=(redundancy>10 + 1),s=40)
        datacursor(display='single',formatter="x={x:.4f}\ny={y:.4f}".format)
        plt.show()


if __name__ == "__main__":
#    HexArray(ditheredInriggers = True)
#    HexArray(AaronsInriggers = False)
#    HexArray(hexNum = 7, AaronsOutriggers = True)
#    HexArray(hexNum = 11, LoadHERAOutriggers = True)
    HexArray(hexNum = 11, AaronsInriggers = True)
#    HexArray(hexNum = 11, AaronsInriggers = True)
    
    