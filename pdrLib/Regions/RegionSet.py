import numpy as np
from Region import Region


class RegionSet:

    def __init__(self, dataMask, image=None):
        self.DataMask = dataMask
        self.Image = image
        self.getParams()
        self.Regions = {}
        for nReg in self.idRegions:
            self.Regions[nReg] = Region(nReg, dataMask, image=image)

    def getParams(self):
        self.idRegions = np.sort(np.unique(self.DataMask))
        if 0 in self.idRegions:
            self.idRegions = self.idRegions[1:]
        self.nRegs = self.idRegions.size

    def getAdjacentRegs(self):
        listOfRegions = self.Regions.keys()
        adjacentRegions = []
        for ii, regionA in enumerate(listOfRegions):
            for regionB in listOfRegions[ii + 1:]:
                if areAdjacent(self.Regions[regionA], self.Regions[regionB]):
                    adjacentRegions.append([regionA, regionB])
        return adjacentRegions

    def getCloseRegs(self, minPixels, adjacent=True):
        closeRegions = []
        listOfRegions = self.Regions.keys()
        for ii, regionA in enumerate(listOfRegions):
            centroidA, radiusA = self.Regions[regionA].getCentroid()
            for regionB in listOfRegions[ii + 1:]:
                centroidB, radiusB = self.Regions[regionB].getCentroid()
                if np.sum((centroidB - centroidA) ** 2) < minPixels ** 2:
                    if adjacent:
                        closeRegions.append([regionA, regionB])
                    else:
                        if not areAdjacent(self.Regions[regionA], self.Regions[regionB]):
                            closeRegions.append([regionA, regionB])
        return closeRegions

    def joinRegions(self, idA, idB):
        self.DataMask[self.Regions[idB].where] = idA
        self.Regions[idA] = Region(idA, self.DataMask, image=self.Image)
        del self.Regions[idB]
        self.getParams()

    def deleteReg(self, id):
        self.DataMask[self.Regions[id].where] = 0
        del self.Regions[id]
        self.getParams()


def areAdjacent(regionA, regionB):
    if regionB.nReg not in regionA.sliceMask:
        return False
    for ii in [-1, 1]:
        for jj in [-1, 1]:
            shiftedArray = (regionA.pixels + np.array([ii, jj])).tolist()
            if True in [item in shiftedArray for item in regionB.pixels.tolist()]:
                return True

    return False
