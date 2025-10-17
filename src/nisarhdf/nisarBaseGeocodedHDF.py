#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 08:44:36 2024

@author: ian
"""


import scipy
import numpy as np
from nisarhdf.nisarBaseHDF import nisarBaseHDF
import pyproj

class nisarBaseGeocodedHDF(nisarBaseHDF):
    '''
    This class creates objects to work with nisar Level 2 geocoded images.
    '''

    def __init__(self,  **keywords):
        '''
        NISAR ROFF HDF Reader

        Parameters
        ----------
        sar : str, optional
            SAR Idendtifier (always LSAR for now). The default is 'LSAR'.
        frequency : str, optional
            frequency band to extract. The default is 'frequencyA'.
        polarization : str, optional
            Polarization. The default is 'HH'.
        layer : str, optional
            Layer to use. The default is 'layer1'.

        Returns
        -------
        None.

        '''
        self.gridType = 'radarGrid'
        self.coordType = 'GEO'
        nisarBaseHDF.__init__(self, **keywords)

        self.levelParams = ['xSize', 'ySize']

    def getEPSG(self):
        '''
        Get epsg
        Returns
        -------
        None.
        '''
        self.epsg = self.toScalar(
            self.h5[self.product]['metadata'][self.gridType]['projection'])
        
    def getPolarizations(self):
        '''
        Parse the list of polarizations
        '''
        self.polarizations = [x.decode("utf-8") for x in 
            self.h5[self.product]['grids'][self.frequency]['listOfPolarizations']]

#
# Interpolation/Data Cube routines
#

    def setupXYGrid(self):
        '''
        Form the X/Y grid for the product

        Returns
        -------
        None.

        '''
        self.xGrid, self.yGrid = np.meshgrid(
            self.xCoords, self.yCoords, indexing='xy')

    def getExtent(self):
        '''
        Get extent (left, right, bottom, top) for imshow

        Returns
        -------
        None.

        '''
        self.extentUL = [self.xCoords[0] - 0.5 * self.dx,
                         self.xCoords[-1] + 0.5 * self.dx,
                         self.yCoords[-1] - 0.5 * self.dy,
                         self.yCoords[0] + 0.5 * self.dy]

        self.extentLL = [self.xCoords[0] - 0.5 * self.dx,
                         self.xCoords[-1] + 0.5 * self.dx,
                         self.yCoords[0] - 0.5 * self.dy,
                         self.yCoords[-1] + 0.5 * self.dy]

    def setupDataCube(self, cubeName, zGridPts=None):
        '''
        Create an interpolator for the incidence angle data cube

        Parameters
        ----------
        cubeName: str
            name of data cube
        Returns
        -------

        interpolator for the cube

        '''
        metadata = self.h5[self.product]['metadata']
        # Adjust for x, y, z order
        cubeData = np.moveaxis(np.moveaxis(
            np.array(metadata['radarGrid'][cubeName]),  0, -1), 0, 1)
        # Get coordinates
        xCoordinates = np.array(metadata['radarGrid']['xCoordinates'])
        yCoordinates = np.array(metadata['radarGrid']['yCoordinates'])
        heightAboveEllipsoid = \
            np.array(metadata['radarGrid']['heightAboveEllipsoid'])
        # Extract subset for baseline
        if zGridPts is not None:
            heightAboveEllipsoid = heightAboveEllipsoid[zGridPts]
        # Create Interpolator
        setattr(self,
                f'{cubeName}Interpolator',
                scipy.interpolate.RegularGridInterpolator(
                        (xCoordinates, yCoordinates, heightAboveEllipsoid),
                        cubeData, bounds_error=False, fill_value=np.nan))

    def interpGrid(self, interpolator, xGrid, yGrid, z, maskNoData=None):
        '''
        Apply regular grid interpolator for x, y coordinates for a date cube.

        Parameters
        ----------
        interpolator : TYPE
            The regularGrid interpolator to use.
        x, y, z : list or nd.array
            x, y, z coordinates to interpolate to.
        maskNoData : bool np.array, optional
            Locations to mask. use instance.noDataLocation when using the full
            product grid (.xGrid and .yGrid). The default is None.
        Returns
        -------
        result.

        '''
        result = interpolator(list(zip(self.ensureNP(xGrid).flatten(),
                                       self.ensureNP(yGrid).flatten(),
                                       self.ensureNP(z).flatten()))
                              ).reshape(self.ensureNP(xGrid).shape)
        if maskNoData is not None:
            # Only compute when needed
            if self.noDataLocations is None and len(self.dataFields) > 0:
                self.noDataLocations = np.isnan(getattr(self, self.dataFields[0]))
            #
            result[self.noDataLocations] = np.nan
        return result

    def rangeDopplerCube(self, x, y, z, maskNoData=None, image='reference',
                         save=False):
        '''
        Interpolate slantRange and zeroDoppler coordinates from data cube

        Parameters
        ----------

        x : nparray
            x coordinates.
        y : nparray
            y coordinates.
        z : nparray
            Elevation coordinates.
        image : str (optional)
            Coordinate are for reference or secondary image.
            Default is 'reference'
        save: bool, optional
            If true, save the result as self.result (e.g.,
            self.slantRangeGrid, self.zeroDopplerTimeGrid). The default is
            False.
        Returns
        -------
        slantRange, zeroDopplerTime : nparrays
             slantRange and zeroDopplerTime coordinates

        '''
        if self.product in ['GCOV']:
            slantRange = 'slantRange'
            zeroDoppler = 'zeroDopplerAzimuthTime'
        else:
            slantRange = f'{image}SlantRange'
            zeroDoppler = f'{image}ZeroDopplerAzimuthTime'
        # Setup interpolators if needed.
        if not hasattr(self, slantRange):
            self.setupDataCube(slantRange, zGridPts=None)
        if not hasattr(self, zeroDoppler):
            self.setupDataCube(zeroDoppler, zGridPts=None)
        # interpolate
        slantRange = self.interpGrid(
            getattr(self, f'{slantRange}Interpolator'),
            x, y, z,
            maskNoData=maskNoData)
        zeroDopplerTime = self.interpGrid(
            getattr(self, f'{zeroDoppler}Interpolator'),
            x, y, z,
            maskNoData=maskNoData)
        self.slantRangeGrid = slantRange
        self.zeroDopplerTimeGrid = zeroDopplerTime
        return slantRange, zeroDopplerTime

#
# I/O routines
#

    def getGeoTransform(self, tiff=True, grimp=False):
        '''
        Compute the geotranform for R/D coordinates

        Returns
        -------
        geotranform : list
            A geotransform with slant range and zero Doppler times

        '''
        xOrigin = self.x0 - self.dx / 2
        
        dy = self.dy
        if grimp and not tiff:
            dy *= -1
            yOrigin = self.y0 - self.dy / 2
        else:
            yOrigin = self.y0 + (self.ySize + 0.5) * np.abs(self.dy)
        #
        geoTransform = [xOrigin, self.dx, 0., yOrigin, 0., dy] 
        if self.downsampleFactorRow > 1  or self.downsampleFactorColumn > 1:
            geoTransform = self.rescale_geoTransform(geoTransform,
                                                     self.downsampleFactorColumn,
                                                     self.downsampleFactorRow)
        return geoTransform
#
# GrIMP Specific Code
#
    def getGeoCoordinates(self):
        '''
        Get Coordinate information (origin, size, and spacing) and asssemble
        information for a geodat file.

        Returns
        -------
        None.

        '''
        # Epsg
        self.epsg = self.toScalar(
            self.h5[self.product]['metadata']['radarGrid']['projection'])
        #
        frequency = self.h5[self.product]['grids'][self.frequency]
        if self.product in ['GUNW', 'GOFF']:
            products = frequency[self.productType][self.polarization]
        elif self.product in ['GCOV']:
            products = frequency
        #
        if self.layer is not None:
            products = products[self.layer]

        # size
        self.xSize = len(products['xCoordinates'])
        self.ySize = len(products['yCoordinates'])
        # Coordinates
        self.xCoords = np.array(products['xCoordinates'])
        self.yCoords = np.array(products['yCoordinates'])
        # Coordinates transform for GrIMP products
        self.x0 = np.min(self.xCoords)
        self.y0 = np.min(self.yCoords)
        self.dx = self.toScalar(products['xCoordinateSpacing'])
        self.dy = self.toScalar(products['yCoordinateSpacing'])
        self.geodat = [[self.xSize, self.ySize],
                       [np.abs(self.dx), np.abs(self.dy)],
                       [self.x0/1000, self.y0/1000.]]
        #
        self.xytoll = pyproj.Transformer.from_crs(f"EPSG:{self.epsg}",
                                                  "EPSG:4326").transform
       
        self.centerLat, self.centerLon = self.xytoll(
            self.x0 + 0.5 * np.abs(self.dx) * self.xSize,
            self.y0 + 0.5 * np.abs(self.dy) * self.ySize)
        #

    def writeGeodat(self, filename):
        '''
        Write geolocation information to a GrIMP geodat file.
        '''
        with open(filename, 'w') as fpOut:
            print('# 2', file=fpOut)
            print(';\n; Image size (pixels) nx ny\n;', file=fpOut)
            print(f'{self.nx:d} {self.ny:d}', file=fpOut)
            print(';\n; Pixel size (m) deltaX deltaY\n;', file=fpOut)
            print(f'{np.abs(self.dx):.4f} {np.abs(self.dy):.4f}', file=fpOut)
            print(';\n; Origin, lower left corner (km) Xo  Yo\n;', file=fpOut)
            print('{self.x0/1000.:.4f} {self.y0/1000.:.4f}', file=fpOut)
            print('&', file=fpOut)
