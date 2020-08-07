from scipy.spatial import distance
import numpy as np
from collections import OrderedDict

class CentroidTrackerAndMitosis():
    def __init__(self):
        # dictionary to track centroid 
        self.centroids = OrderedDict()
        # dictionary to track coordiantes
        self.allCoordiantes = OrderedDict()
        # value of new unique id for cell 
        self.nextCellID = 0
        # dictionary to track lost cells in case of bounding box glitch for only one frame
        self.lost = OrderedDict()
        # area of cell
        self.area = OrderedDict()
        self.mitosis = []

    def register(self, centroid, coordinates, area):
        # register centroid of cell with next available id
        self.centroids[self.nextCellID] = centroid
        # register coordinate of cell with next available id
        self.allCoordiantes[self.nextCellID] = coordinates
        # no lost 
        self.lost[self.nextCellID] = 0
        # register area of cell 
        self.area[self.nextCellID] = area
        # increment cell id 
        self.nextCellID += 1
        
    def remove(self, cellId):
        # cell no longer exists 
        # delete centroid value of that cell in dictionary
        del self.centroids[cellId]
        # delete coordinate values of that cell in dictionary
        del self.allCoordiantes[cellId]
        # delete lost value of that cell in dictionary
        del self.lost[cellId]
        # delete area value of that cell in dictionary 
        del self.area[cellId]    

    def update(self, parameters, dataset):
    
        # array to store current centroid vals
        currentCentroids = np.zeros((len(parameters), 2), dtype="int")
        
        # array to store current coordinate vals from parameters
        currentCoordiantes = np.zeros((len(parameters), 4), dtype="int")
        
        # array to store current area vals of cell
        currentArea = np.zeros((len(parameters), 1), dtype="int")
  
        # get value of centroid and coordinates 
        for (i, (wighted_center_x, wighted_center_y, area, X_0, Y_0, X_n, Y_n)) in enumerate(parameters):
            currentCentroids[i] = (int(wighted_center_x),int(wighted_center_y))
            currentCoordiantes[i] = (X_0, Y_0,X_n, Y_n)
            currentArea[i] = area
            
        # for case where no cells are being tracked 
        if len(self.centroids) is 0:
            for i in range(0, len(currentCentroids)):
                # register centroids and (x,y)
                self.register(currentCentroids[i], currentCoordiantes[i], currentArea[i] )

        # for case when we need to match new centroids with existing centroids 
        else:
            # existing centroids 
            cellCentroids = np.array(list(self.centroids.values()))
            
            #existing ids 
            cellIds = list(self.centroids.keys())

            # compute each distance 
            # totalDistance is the distance matrix of all centroid distance
            totalDistance = distance.cdist(cellCentroids, currentCentroids, metric='euclidean')

            # sort by minimum row and then use that to sort collumns 
            rows = totalDistance.min(axis=1).argsort()
            cols = totalDistance.argmin(axis=1)[rows]

            # track of used columns and rows 
            rowsHistory, colsHistory = set(), set()
            fatherMit = self.mitosis
            self.mitosis = []
            fatherCentroids = self.centroids
            for (row, col) in zip(rows, cols):
                # already checked 
                if row in rowsHistory or col in colsHistory:
                    continue

                # get id of current row and update its new centroid 
                # reset lost value
                cellId = cellIds[row]
                self.centroids[cellId] = currentCentroids[col]
                self.allCoordiantes[cellId] = currentCoordiantes[col]
                #self.aread[cellId] = currentArea[col]
                #(X_0, Y_0, X_n, Y_n) = currentCoordiantes[col]
                oldArea = self.area[cellId]
                self.area[cellId] = currentArea[col]
                self.lost[cellId] = 0
        
                # add these row and col to history
                colsHistory.add(col)
                rowsHistory.add(row)
                if oldArea > 0 and self.area[cellId] > 0:
                    ratio = oldArea/self.area[cellId] 
                else:
                    ratio = 1
                if dataset == 2:
                    if (ratio > 1.6 or ratio < 0.4 ):
                        print("Predicted Mitosis " + str(cellId) + ":" + str(self.centroids[cellId]))
                        # found possible mitosis
                        (x,y) = self.centroids[cellId]
                        self.mitosis.append((cellId,x,y))
                elif dataset == 3:
                    if (ratio > 1.5 or ratio < 0.5 ):
                        print("Predicted Mitosis " + str(cellId) + ":" + str(self.centroids[cellId]))
                        # found possible mitosis
                        (x,y) = self.centroids[cellId]
                        self.mitosis.append((cellId,x,y))
                elif dataset == 1:
                    if (ratio > 1.7 or ratio < 0.3 ):
                        print("Predicted Mitosis " + str(cellId) + ":" + str(self.centroids[cellId]))
                        # found possible mitosis
                        (x,y) = self.centroids[cellId]
                        self.mitosis.append((cellId,x,y))
                else:
                    print("dataset has not been entered")
                    if (ratio > 1.6 or ratio < 0.4 ):
                        print("Predicted Mitosis " + str(cellId) + ":" + str(self.centroids[cellId]))
                        # found possible mitosis
                        (x,y) = self.centroids[cellId]
                        self.mitosis.append((cellId,x,y))
            print("")
            # for the rows and cols we havnt done anything to 
            rowsToCheck = set(range(0, totalDistance.shape[0])).difference(rowsHistory)
            colsToCheck = set(range(0, totalDistance.shape[1])).difference(colsHistory)
            
            
            # if number of cell centroid is greater than number of new centroid 
            # check if cells are lost 
            if totalDistance.shape[0] > totalDistance.shape[1]:
                for row in rowsToCheck:
                    cellId = cellIds[row]
                    self.lost[cellId] += 1

                    if self.lost[cellId] >= 1:
                        self.remove(cellId)

            # register any new centroids and coordinates 
            else:
                for col in colsToCheck:
                    self.register(currentCentroids[col], currentCoordiantes[col], currentArea[col])
      
            if len(fatherMit) != 0:
                for (mcellId,x,y) in fatherMit:
                    if mcellId in self.centroids.keys():
                        theCellInMitosis = np.array([[x,y]])
                        cellCentroidsForMitosis = np.array(list(self.centroids.values()))
                        mitosisDistance = distance.cdist(theCellInMitosis, cellCentroidsForMitosis, metric='euclidean')
                        sortedMitosisDistance = sorted(mitosisDistance[0])
                        index=np.where(mitosisDistance[0] == sortedMitosisDistance[1])
                        indexValue = index[0]
                        sibilingCellIdValue = cellCentroidsForMitosis[indexValue[0]]

                        print("Mother Cell's ID is "+ str(mcellId) + ":" + str((x,y)))
                        self.register(self.centroids[mcellId], self.allCoordiantes[mcellId], self.area[mcellId])
                        self.remove(mcellId)
                        theID = None
                        for (i,j) in self.centroids.items():
                            if (j == sibilingCellIdValue).all():
                                theID = i
                                print("Daughter Cell 1 ID is " + str(theID) + ":" + str(self.centroids[theID]))
                                print("Daughter Cell 2 ID is " + str(self.nextCellID-1) + ":" + str(self.centroids[self.nextCellID-1])+"\n")
                                for (IDToComp,x1,y1) in self.mitosis:
                                    if mcellId == IDToComp:
                                        self.mitosis.remove((mcellId,x1,y1))
                                        self.mitosis.append((self.nextCellID-1,x1,y1))  
                    
        # return centroids and coordinates 
        return (self.centroids, self.allCoordiantes, self.area, self.mitosis)