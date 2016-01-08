# GENERATOR OF HEXAGONAL ARRAYS WITH INRIGGERS AND OUTRIGGERS FOR ARRAY COFINGURATION PROJECT
# by Josh Dillon
#
# This code generates a file of antenna positions and a files with baselines, baseline info, and redundancies

import sys,os,os.path
import numpy as np
import cPickle as pickle
import matplotlib
import matplotlib.pyplot as plt
#plt.close('all')

def HexArray(Separation = 14.6, hexNum = 11, 
             SplitCore = False,
             NineWaySplitCore = False,
             SplitCoreOutriggers = False,
             JoshsOutriggers = False, 
             RedundantOutriggers = False, 
             LoadHERAOutriggers = False,
             fiducialInriggers = False, 
             redundantTriangleInriggers = False,
             redundantPairInriggers = False,
             redundantHexInriggers = False,
             halfCornerInriggers = False,
             fullCornerInriggers = False, **kwargs):

    precisionFactor = 1000000        
    
    #Calculating Positions:
    #Main Hex
    positions = [];
#    for row in range(hexNum-1,-(hexNum)+SplitCore,-1):
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

#    #Split the core into 3 pieces
    if SplitCore:
        for i,pos in enumerate(positions):          
            theta = np.arctan2(pos[1],pos[0])
            if not (pos[0]==0 and pos[1]==0):
                if (theta > -np.pi/3 and theta < np.pi/3):
                    positions[i] = np.asarray(pos) + (upRight + upLeft)/3
                if (theta >= np.pi/3 and theta < np.pi):
                    positions[i] = np.asarray(pos) +upLeft  - (upRight + upLeft)/3



    if NineWaySplitCore:
        from sympy.geometry import RegularPolygon, Point, Polygon
        distances = np.asarray([np.linalg.norm(pos) for pos in positions])
        fullRadius = np.max(distances)
        interiorHex = RegularPolygon(Point(0, 0), 1.001*3**-.5*fullRadius, 6)
        interiorHex.spin(np.pi/6)
        interiorHex = Polygon(*[Point(vert.x.evalf(), vert.y.evalf()) for vert in interiorHex.vertices])

        positions.remove(positions[76])
    
        isInterior = np.zeros(len(positions),dtype=bool)    
        for i,pos in enumerate(positions):
            if distances[i] < fullRadius/2.0: isInterior[i] = True
            elif distances[i] <= 3**-.5*fullRadius and interiorHex.encloses_point(Point(pos[0],pos[1])): isInterior[i] = True
        



        v1 = (upRight + upLeft)/3       
        v2 = right/3
        
        for i,pos in enumerate(positions):
            theta = np.arctan2(pos[1],pos[0])
            if isInterior[i]:
                if not (pos[0]==0 and pos[1]==0):
                    if (theta > -np.pi/3 and theta < np.pi/3): #right
                        positions[i] = np.asarray(pos) +v1
                    if (theta >= np.pi/3 and theta < np.pi): #top left
                        positions[i] = np.asarray(pos) +upLeft  - v1
            else:
                if (theta > -np.pi/6 and theta < np.pi/6): #right
                    positions[i] = np.asarray(pos) + 2*v2
                if (theta <= -5*np.pi/6 or theta >= 5*np.pi/6): #left
                    positions[i] = np.asarray(pos) -v2 - right
                if (theta > -np.pi/2 and theta <= -np.pi/6): #bottom right
                    positions[i] = np.asarray(pos) -v1 +v2
                if (theta > -5*np.pi/6 and theta <= -np.pi/2): #bottom left
                    positions[i] = np.asarray(pos) -v1-v2
                if (theta >= np.pi/6 and theta < np.pi/2): #top right
                    positions[i] = np.asarray(pos) + v1 + 2*v2
                if (theta >= np.pi/2 and theta < 5*np.pi/6): #top left
                    positions[i] = np.asarray(pos) +v1-2*v2 

    nCore = len(positions)

    #Aaron's fiducial inriggers
    if fiducialInriggers:
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
            inriggerCenters = [((hexNum-1)/2) * (upLeft + upRight) + upLeft + upLeft + right + (right + upRight)/3,       #top
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
     
    #Redundantly calibratable inriggers with full UV coverage
    if fullCornerInriggers or halfCornerInriggers:
        positions.append(hexNum * (upLeft) + (right + upRight)/3)
        positions.append(hexNum * (upLeft) - (right + upRight)/3)
        positions.append(hexNum * (upRight) + (-right + upLeft)/3)
        positions.append(hexNum * (upRight) - (-right + upLeft)/3)
        positions.append(hexNum * (right) + (upRight + upLeft)/3)
        positions.append(hexNum * (right) - (upRight + upLeft)/3)        
        if fullCornerInriggers:
            positions.append(-(hexNum * (upLeft) + (right + upRight)/3))
            positions.append(-(hexNum * (upLeft)  - (right + upRight)/3))
            positions.append(-(hexNum * (upRight) + (-right + upLeft)/3))
            positions.append(-(hexNum * (upRight) - (-right + upLeft)/3))
            positions.append(-(hexNum * (right) + (upRight + upLeft)/3))
            positions.append(-(hexNum * (right) - (upRight + upLeft)/3))



    #Half Wavelength Stuff
#    for i in range(1,14):
#        positions.append(i*(upLeft + upRight)/27)
        
 


#%%    
    
    #Outriggers
    
    if SplitCoreOutriggers:
        exteriorHexNum = 4
        for row in range(exteriorHexNum-1,-(exteriorHexNum),-1):
            for col in range(2*exteriorHexNum-abs(row)-1):
                xPos = ((-(2*exteriorHexNum-abs(row))+2)/2.0 + col)*Separation*(hexNum-1)
                yPos = row*Separation*(hexNum-1)*3**.5/2
                theta = np.arctan2(yPos,xPos)       
                if ((xPos**2 + yPos**2)**.5 > Separation*(hexNum+1)):
                    if (theta > -np.pi/3 and theta < np.pi/3):
                        positions.append(np.asarray([xPos, yPos, 0]) + (upRight + upLeft)/3)
                    elif (theta >= np.pi/3 and theta < np.pi):
                        positions.append(np.asarray([xPos, yPos, 0]) +upLeft - (upRight + upLeft)/3)
                    else:                    
                        positions.append([xPos, yPos, 0])


    if JoshsOutriggers:
        outriggerHexNum = 3
        for row in range(outriggerHexNum-1,-(outriggerHexNum),-1):
            for col in range(2*outriggerHexNum-abs(row)-1):
                xPos = ((-(2*outriggerHexNum-abs(row))+2)/2.0 + col)*Separation*(np.floor(1.5*hexNum));
                yPos = row*Separation*(np.floor(1.5*hexNum))*3**.5/2;
                if xPos != 0 or yPos != 0:    
                    positions.append([xPos, yPos, 0])

    if RedundantOutriggers:
        outriggerHexNum = 3
        for row in range(outriggerHexNum-1,-(outriggerHexNum),-1):
            for col in range(2*outriggerHexNum-abs(row)-1):
                yPos = ((-(2*outriggerHexNum-abs(row))+2)/2.0 + col)*Separation*(3**.5*(hexNum-1));
                xPos = row*Separation*(3**.5*(hexNum-1))*3**.5/2;
                if xPos != 0 or yPos != 0:    
                    positions.append([xPos, yPos, 0])
        
    if LoadHERAOutriggers:
        outriggerPositions = np.loadtxt(os.path.dirname(os.path.abspath(__file__)) + '/hexconfig352_outriggers_only.dat')
        for outrigger in outriggerPositions:
            positions.append(np.append(outrigger,0))
        precisionFactor = 10000

    #Calculating Baselines    
    nAntennas = len(positions)
    baselines = []
    baselinePairs = []
    #print "WARNING: DOUBLING UNIQUE BASELINES!!!"
    for ant1 in range(nAntennas):
        for ant2 in range(ant1+1,nAntennas):
            deltax = int(np.round(precisionFactor*(positions[ant1][0]-positions[ant2][0])))
            deltay = int(np.round(precisionFactor*(positions[ant1][1]-positions[ant2][1])))
            deltaz = int(np.round(precisionFactor*(positions[ant1][2]-positions[ant2][2])))
            if (deltax**2+deltay**2+deltaz**2 < precisionFactor**2*Separation**2): print "WARNING: antennas " + str(ant1) + " and " + str(ant2) + " are too close together!"
            if deltay > 0 or (deltay == 0 and deltax > 0):                
                baselines.append((deltax, deltay, deltaz))
                baselinePairs.append((ant1, ant2))
            else:
                baselines.append((-deltax, -deltay, -deltaz))
                baselinePairs.append((ant2, ant1))                

    print len(baselinePairs)    
    
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
        from mpldatacursor import datacursor        
        uniqueBaselines = np.asarray([uniqueBaseline[0] for uniqueBaseline in baselineDict.items()])/(1.0*precisionFactor)
        redundancy = np.asarray([len(uniqueBaseline[1]) for uniqueBaseline in baselineDict.items()])
        uniqueBaselines = np.append(uniqueBaselines, -uniqueBaselines, axis=0)
        redundancy = np.append(redundancy, redundancy, axis=0)
        plt.figure(1)
        plt.clf()
        plt.scatter(np.asarray(positions)[:,0]/Separation,np.asarray(positions)[:,1]/Separation)
        plt.axis('equal')
        plt.figure(2)    
        plt.clf()
#        plt.scatter(uniqueBaselines[:,0]/1.0/Separation, uniqueBaselines[:,1]/1.0/Separation,c=np.minimum(redundancy,100000),s=40)
        plt.scatter(uniqueBaselines[:,0]/1.0/1.5, uniqueBaselines[:,1]/1.0/1.5,c=np.minimum(redundancy,100000),s=40)
        plt.colorbar()
        plt.title('Baseline Redundancy')
        #datacursor(display='single',formatter="x={x:.4f}\ny={y:.4f}".format)
        plt.show()
        plt.axis('equal')
        return np.asarray(positions)

if __name__ == "__main__":
 #   positions = HexArray(hexNum = 7, RedundantOutriggers = False, fullCornerInriggers = True)
#    positions = HexArray(hexNum = 7, SplitCore = False, SplitCoreOutriggers = True)
#    positions = HexArray(hexNum = 11, NineWaySplitCore = True, SplitCoreOutriggers = True)
    positions = HexArray(hexNum = 7)

    
    
    if False:
        import omnical.info as oi
        import omnical.arrayinfo as oai
        redundantInfo = oi.RedundantInfo()
        reds = oai.compute_reds(positions)
        redundantInfo.init_from_reds(reds, positions)
        AtransA = redundantInfo.At.dot(redundantInfo.At.T).toarray()
        BtransB = redundantInfo.Bt.dot(redundantInfo.Bt.T).toarray()    
        gainVariances = np.diag(np.linalg.pinv(AtransA)[0:len(positions),0:len(positions)])
        gainVariances = gainVariances/gainVariances.min()
        print "Gain and phase modes that can't be Omnicaled:"
#        print [np.sum(np.abs(np.linalg.eigvals(XtransX)<1e-10)) for XtransX in [AtransA, BtransB]]
        print [len(XtransX) - np.linalg.matrix_rank(XtransX, tol=1e-10) for XtransX in [AtransA, BtransB]]
        plt.figure(3)
        plt.clf()
        hexScatter = plt.scatter(positions[:,0], positions[:,1], c=np.sqrt(gainVariances), s=100)
        plt.colorbar()
        plt.title('Antenna Relative Gain Errors')
