import numpy as np

class Region:
    
    def __init__(self, nReg, dataMask, image=None):
        self.nReg = nReg
        self.where = np.where(dataMask == nReg)
        self.pixels = np.array([self.where[0], self.where[1]]).T
        self.dataMask = dataMask
        self.sliceMask, self.sliceMaskWhere = self.getSlice(self.dataMask)
        # self.contour = self.getContour()
        # self.nearbyRegs = self.getNearbyRegs()
        if image != None:
            self.sliceImage, dummy = self.getSlice(image)
            self.image = image 
            self.peak = np.max(image[self.where])
            peakWhere = image[self.where].argmax()
            self.peakWhere = (np.array(self.where[0][peakWhere]), 
                np.array(self.where[1][peakWhere]))
    
    def getSlice(self, dataMask, extend=1):
        dims = dataMask.shape
        iMin, iMax = self.pixels[[0,-1], 0]
        jMin, jMax = np.sort(self.pixels[:,1])[[0,-1]]
        try:
            slice = dataMask[iMin-extend:iMax+(1+extend), jMin-extend:jMax+(1+extend)]
        except IndexError as detail:
            iSize = iMax-iMin+1
            jSize = jMax-jMin+1
            slice = np.zeros([iSize+2*extend, jSize+2*extend], dtype=int)
            slice[extend:iSize+extend, extend:jSize+extend] = dataMask[iMin:iMax+1, jMin:jMax+1]
            if iMin > extend-1: 
                slice[0, 1:jSize+1] = dataMask[iMin-extend:iMin, jMin:jMax+1]
            if iMax < dims[0]-extend: 
                slice[iSize+1, 1:jSize+1] = dataMask[iMax+1:iMax+extend, jMin:jMax+1]
            if jMin > extend-1: 
                slice[1:iSize+1, 0] = dataMask[iMin:iMax+1, jMin-extend:jMin]
            if jMax < dims[1]-extend: 
                slice[1:iSize+1, jSize+1] = dataMask[iMin:iMax+1, jMax+1:jMax+extend]
            # if iMin > 0 and jMin > 0: 
            #     slice[0, 0] = dataMask[iMin-1, jMin-1]
            # if iMin > 0 and jMax < dims[1]-1 : 
            #     slice[0, jSize+1] = dataMask[iMin-1, jMax+1]
            # if iMax < dims[0]-1 and jMin > 0: 
            #     slice[iSize+1, 0] = dataMask[iMax+1, jMin-1]
            # if iMax < dims[0]-1 and jMax < dims[1]-1: 
            #     slice[iSize+1, jSize+1] = dataMask[iMax+1, jMax+1]
        return slice, np.where(slice==self.nReg)
    
    
    def getContour(self):
        diff = np.array([self.pixels[0,0], np.sort(self.pixels[:,1])[0]]) - 1
        pix = np.where(self.sliceReg == self.nReg)
        sum = 0
        for ii in [-1,0,1]:
            for jj in [-1,0,1]:
                sum += self.sliceReg[pix[0]+ii, pix[1]+jj]
        noBorder = np.where(sum/9. != self.nReg)
        return diff[0] + pix[0][noBorder], diff[1] + pix[1][noBorder]
    
    
    def getNearbyRegs(self):
        regions = np.unique(self.sliceMask)
        regions = np.delete(regions, np.where((regions==0) | (regions==self.nReg)))
        return regions
    
    
    # def addData(self, image):
        # self.sliceData, dummy = self.getSlice(image)
        # return
    
    
    def getFlux(self, convFactor=1., image=None):
        if image == None:
            if hasattr(self, 'sliceImage'):
                dataTmp = self.sliceImage[self.sliceMaskWhere]
            else:
                raise NameError('Image data not available')
        else:
            dataTmp = image[self.where]
        self.flux = np.sum(dataTmp) * convFactor
        return self.flux
    
    
    # def getMax(self, convFactor=1., data=None):
    #     if data == None:
    #         if hasattr(self, 'sliceData'):
    #             dataTmp = self.sliceData[self.sliceMaskWhere]
    #         else:
    #             raise NameError('Image data not available')
    #     else:
    #         dataTmp = data[self.where]
    #     self.max = np.max(dataTmp) * convFactor
    #     return self.max
    
    
    def getCentroid(self):
        centroid = np.array([np.mean(self.pixels[:,0]), np.mean(self.pixels[:,1])])
        radius = np.sqrt(np.max(np.sum((self.pixels-centroid)**2, axis=1)))
        self.centroid = centroid
        self.radius = radius
        return centroid, radius
    
    