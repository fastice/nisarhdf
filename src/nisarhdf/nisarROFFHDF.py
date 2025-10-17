#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 08:36:38 2024

@author: ian
"""
from nisarhdf.nisarBaseRangeDopplerHDF import nisarBaseRangeDopplerHDF
from osgeo import gdal
import os
import numpy as np
import rioxarray
from nisarhdf import writeMultiBandVrt


class nisarROFFHDF(nisarBaseRangeDopplerHDF):
    '''
    This class creates objects to work with nisar RUNWHDF images.
    '''

    def __init__(self, sar='LSAR',
                 frequency='frequencyA', polarization=None,
                 layer='layer1', byteOrder='LSB', productType='pixelOffsets'):
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
        nisarBaseRangeDopplerHDF.__init__(self,
                                          sar=sar,
                                          product='ROFF',
                                          frequency=frequency,
                                          productType=productType,
                                          polarization=None,
                                          layer=layer,
                                          productData='alongTrackOffset',
                                          bands='swaths',
                                          bytOrder=byteOrder)
        self.lookType = 'Offset'
        self.isSecondary = False
        self.productParams = ['r0', 'a0', 'deltaR', 'deltaA']
        for param in self.RDParams:
            self.productParams.append(f'{self.lookType}{param}')

    def parseParams(self, productType='pixelOffsets', polarization=None,
                    layers=['layer1', 'layer2', 'layer3'],
                    fields=None, secondGeodat=None, noLoadData=False,
                    **keywords):
        '''
        Parse all the params needed for offsets

        Returns
        -------
        None.

        '''
        self.getPolarization(polarization)
        self.getOrbitAndFrame(**keywords)
        self.getLookDirection()
        self.getRangeBandWidth()
        self.getOrbitPassDirection()
        self.parseRefDate()
        self.getOffsetSize()
        self.getOffsetWindowParams()
        self.getOffsetParams()
        self.getSLCSlantRange()
        self.getSLCZeroDopplerTime()
        self.getOffsetSlantRange()
        self.getOffsetZeroDopplerTime()
        self.effectivePRF()
        self.getCenterIncidenceAngle()
        self.getEPSG()
        self.getCorners()
        self.getCenterLatLon()
        # print(self.corners)
        self.getWavelength()
        self.getGranuleNames()
        self.SLCSceneCenterAlongTrackSpacing()
        self.getExtent()
        self.fieldDict = {'pixelOffsets': ['slantRangeOffset', 
                                           'slantRangeOffsetVariance',
                                           'alongTrackOffset',
                                           'alongTrackOffsetVariance',
                                           'crossOffsetVariance',
                                           'correlationSurfacePeak',
                                           'snr']}
        #
        # This mess is to handle both layers and things that are not layers
        # So far this is just the DEM
        singleFields = []
        if fields is None:
            layerFields = self.fieldDict['pixelOffsets'] 
            fields = layerFields + ['digitalElevationModel']
            singleFields = ['digitalElevationModel']
        else:
            layerFields = [x for x in fields if x not in ['digitalElevationModel']]
            singleFields = [x for x in fields if x in ['digitalElevationModel']]
        
        # Load data
        self.getLayers(layerFields, layers=layers, noLoadData=noLoadData)
        self.loadData(singleFields, noLoadData=noLoadData,
                      resetFields=False)
        # scale factors to convert back to pixels
        self.scaleFactors = {'slantRangeOffset': 1./self.SLCRangePixelSize,
                             'alongTrackOffset': 1./self.SLCAzimuthPixelSize,
                             'slantRangeOffsetVariance':
                                 (1./self.SLCRangePixelSize)**2,
                             'alongTrackOffsetVariance':
                                 (1./self.SLCAzimuthPixelSize)**2}
        # Suffixes for outputs
        self.suffixes = {'slantRangeOffset': '.dr',
                         'slantRangeOffsetVariance': '.vr',
                         'alongTrackOffset': '.da',
                         'alongTrackOffsetVariance': '.va',
                         'crossOffsetVariance': '.vc',
                         'correlationSurfacePeak': '.cc',
                         'snr': '.snr'}

    def getOffsetZeroDopplerTime(self, secondary=False):
        '''
        Input y5 and return dictionary with zero Doppler spacing and
        first, middle, and last zero Doppler.

        Returns
        -------
        None.

        '''
        self.OffsetZeroDopplerTimeSpacing = \
            self.SLCZeroDopplerTimeSpacing * self.deltaA
        self.getSLCZeroDopplerTime(secondary=secondary)
        # Compute first, middle and last Doppler times
        self.OffsetFirstZeroDopplerTime = self.SLCFirstZeroDopplerTime + \
            self.a0 * self.SLCZeroDopplerTimeSpacing
        # last
        self.OffsetLastZeroDopplerTime = self.OffsetFirstZeroDopplerTime + \
            self.OffsetZeroDopplerTimeSpacing * (self.OffsetAzimuthSize-1)
        # Middle
        self.OffsetMidZeroDopplerTime = 0.5 * \
            (self.OffsetFirstZeroDopplerTime + self.OffsetLastZeroDopplerTime)
        bands = self.h5[self.product][self.bands]
        self.zeroDopplerTimeData = np.array(
            bands[self.frequency][self.productType]['zeroDopplerTime'])
        self.OffsetAzimuthPixelSize = self.SLCAzimuthPixelSize * self.deltaA

    def getOffsetSlantRange(self):
        '''
        Get slant ranges for offsets

        Returns
        -------
        None.

        '''
        self.OffsetRangePixelSize = self.SLCRangePixelSize * self.deltaR
        self.OffsetNearRange = self.SLCNearRange + \
            self.SLCRangePixelSize * self.r0
        self.OffsetFarRange = self.OffsetNearRange + \
            self.OffsetRangePixelSize * (self.OffsetRangeSize - 1)
        self.OffsetCenterRange = 0.5 * \
            (self.OffsetNearRange + self.OffsetFarRange)
        #
        bands = self.h5[self.product][self.bands]
        self.slantRangeData = np.array(
            bands[self.frequency][self.productType]['slantRange'])

    def getOffsetSize(self):
        '''
        Get size of offset products
        '''
        frequency = self.h5[self.product][self.bands][self.frequency]
        polarization = frequency[self.productType][self.polarization]
        self.OffsetAzimuthSize, self.OffsetRangeSize = \
            polarization['layer1']['alongTrackOffset'].shape

    def getOffsetParams(self):
        '''
        Get offset params, a0, r0, deltaA, deltaR

        Returns
        -------
        None.

        '''
        procInfo = self.h5[self.product]['metadata']['processingInformation']
        frequency = procInfo['parameters'][self.productType][self.frequency]
        #
        self.deltaA = np.array(frequency['alongTrackSkipWindowSize']).item()
        self.deltaR = np.array(frequency['slantRangeSkipWindowSize']).item()
        self.a0 = np.array(frequency['alongTrackStartPixel']).item()
        self.r0 = np.array(frequency['slantRangeStartPixel']).item()

    def writeOffsetsVrt(self, vrtFile, sourceFiles, geodat1, geodat2, mask,
                        Image1=None, Image2=None, scaleFactor=1,
                        descriptions=None, byteOrder='LSB',
                        radarCoordinates=False):
        '''
        Write a VRT file to accompany offset products

        Parameters
        ----------
        vrtFile : str
            Name of vrt file.
        sourceFiles : list of str
            List of filenames.
        geodat1 : str
            Name of first geodat to include in meta data.
        geodat2 : str
            Name of second geodat to include in meta data.
        mask : str
            Name of mask.
        Image1 : str, optional
            File name for first image in pair. The default is None.
        Image2 : TYPE, optional
            File name for second image in pair. The default is None.
        scaleFactor : float, optional
            Scale factor to apply to image. The default is 1.
        descriptions : list of str, optional
            Descriptiosn for each band. The default is None.
        byteOrder : str, optional
            Byte order. The default is 'LSB'.
        radarCoordinates : bool, optional
            Set True to use radar coordinates for the geo transform. The
            default is False.

        Returns
        -------
        None.

        '''
        metaData = {}
        #
        if descriptions is None:
            descriptions = [os.path.basename(x) for x in sourceFiles]
        #
        for key in ['a0', 'r0', 'deltaA', 'deltaR']:
            metaData[key] = getattr(self, key)

        for key, var in zip(['geo1', 'geo2', 'mask', 'Image1', 'Image2'],
                            [geodat1, geodat2, mask, Image1, Image2]):
            if var is not None:
                metaData[key] = var
        #
        if radarCoordinates:
            geoTransform = self.rangeDopplerGeoTransform()
        else:
            geoTransform = [-0.5, 1., 0., -0.5, 0., 1.]
        #
        writeMultiBandVrt(vrtFile,
                          self.OffsetRangeSize,
                          self.OffsetAzimuthSize,
                          sourceFiles,
                          descriptions,
                          byteOrder=byteOrder,
                          eType=gdal.GDT_Float32,
                          geoTransform=geoTransform,
                          metaData=metaData,
                          noDataValue=-2.e9)

    def writeOffsetsDatFile(self, datFile, geodat1=None, geodat2=None):
        '''
        Write an offsets dat file (for GrIMP processing only)

        Parameters
        ----------
        datFile : str
            File name for data file.
        geodat1 : str
            Name of first geodat to include in meta data.
        geodat2 : str
            Name of second geodat to include in meta data.

        Returns
        -------
        None.

        '''
        with open(datFile, 'w') as fp:
            print(self.a0, self.r0,
                  self.OffsetRangeSize, self.OffsetAzimuthSize,
                  self.deltaR, self.deltaA,
                  file=fp)
            if geodat1 is None:
                return
            if geodat2 is None:
                geodat2 = geodat1.replace('.in', '.secondary'
                                          ).replace('.geojson', '.secondary')
            print(geodat1, geodat2, file=fp)

    def applyMask(self, maskFile, layers=[3], maskValue=1):
        '''
        Apply mask contained in file

        Parameters
        ----------
        maskFile : str
            vrt file with mask (offsets.mask.vrt)
        t
        layers : list, optional
            List of layers to process. The default is [3].

        Returns
        -------
        None.

        '''
        #
        mask = rioxarray.open_rasterio(maskFile)
        masked = np.squeeze(mask.data) > maskValue
        # Load the product fields if not already done
        self.getProducts(['alongTrackOffset', 'slantRangeOffset'])
        # Loop over layers
        for layer in layers:
            self.alongTrackOffset[layer-1][masked] = np.nan
            self.slantRangeOffset[layer-1][masked] = np.nan

    def rangeDopplerGrid(self):
        '''
        Compute range (m) zero doppler time (sec) coordinates for the offsets.
        grid.

        Returns
        -------
        None.

        '''
        pass

