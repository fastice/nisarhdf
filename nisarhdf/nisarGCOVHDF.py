#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 08:51:41 2024

@author: ian
"""
import os
import h5py
import numpy as np
from nisarhdf import nisarBaseGeocodedHDF


class nisarGCOVHDF(nisarBaseGeocodedHDF):
    '''
    This class creates objects to work with nisar RUNWHDF images.
    '''

    def __init__(self,  sar='LSAR', product='GCOV', frequency='frequencyA',
                 productType=None, polarization=None,
                 layer=None, productData=None, bands='grids',
                 referenceOrbitXML=None, secondaryOrbitXML=None,
                 isSecondary=False):
        '''
        Instantiate nisarVel object. Possible bands are 'image', 'sigma0',
        'gamma0', or user defined.
        Parameters
        ----------
        imageType: str
            imageType custom name, or image, sigma0, gamma0. If not specified
            determined from grimp product name.
        verbose : bool, optional
            Increase level of informational messages. The default is True.
        noData : scalar
            no data value. Defaults to np.nan if not image/sigma/gamma
        Returns
        -------
        None.
        '''
        nisarBaseGeocodedHDF.__init__(self,
                                      sar=sar,
                                      product=product,
                                      frequency=frequency,
                                      productType=productType,
                                      polarization=polarization,
                                      layer=layer,
                                      productData=productData,
                                      bands=bands,
                                      isSecondary=isSecondary,
                                      referenceOrbitXML=referenceOrbitXML,
                                      secondaryOrbitXML=secondaryOrbitXML)
        self.productType = None
        self.productParams = ['NumberRangeLooks', 'NumberAzimuthLooks']
        self.lookType = None

    def parseParams(self, secondary=False, noLoadData=False, sigma0=False,
                    dB=False, fields=None, productType=None, factorFrom=None,
                    **keywords):
        '''
        Parse all the params needed to make a geodatNRxNA.geojson file

        Returns
        -------
        None.

        '''
        #
        self.polarization = None
        self.getOrbitAndFrame(**keywords)     
        self.getNumberOfLooks()
        self.getRangeBandWidth()
        self.getLookDirection()
        self.getOrbitPassDirection()
        self.parseRefDate()
        self.getGeoCoordinates()
        self.getWavelength()
        self.getGranuleNames()
        self.getEPSG()
        self.getSLCSlantRange()
        self.getSLCZeroDopplerTime()
        self.effectivePRF()
        self.getExtent()
        self.covTerms = [self.parseString(x) for x in
                         self.h5[self.product][self.bands][self.frequency][
                        'listOfCovarianceTerms']]
        self.units = 'power'
        #
        # Default to all fields the granule holds. A slim GCOV
        # (downloadNISARoptimized) keeps one covariance term, drops
        # numberOfLooks, and stores rtcGammaToSigmaFactor in a shared file
        # (read with factorFrom), but still lists every term.
        grid = self.h5[self.product][self.bands][self.frequency]
        if fields is None:
            fields = [x for x in self.covTerms +
                      ['mask', 'numberOfLooks', 'rtcGammaToSigmaFactor']
                      if x in grid]
            if factorFrom is not None and 'rtcGammaToSigmaFactor' not in grid:
                fields.append('rtcGammaToSigmaFactor')
        self.factorFrom = factorFrom
        #
        # Remove polarizations that were requested but not available 
        new_fields = []
        missing = []
        for field in fields:
            if field in ['HHHH', 'VVVV', 'HVHV', 'VHVH'] and \
                    (field not in self.covTerms or field not in grid):
               missing.append(field)
            else:
                new_fields.append(field)
        #if len(missing) > 0:
        #     print(f'\033[1mPolarization(s) not in current file\033[0m\n {missing}')
        fields = new_fields
        #
        self.fields = fields
        self.dB = dB
        self.backscatterType = 'gamma0'
        self.loadData(fields, noLoadData=noLoadData)
        #
        self.dB = dB
        if sigma0 and not noLoadData:
            if getattr(self, 'rtcGammaToSigmaFactor', None) is None:
                self.getImageData('rtcGammaToSigmaFactor', useNumpy=True)
            print('computing sigma0')
            self.computeSigma()
        if self.dB and not noLoadData:
            print('computing DB')
            self.computedB()

    def getImageData(self, productField, layer=None, useNumpy=True,
                     power=True):
        '''
        As nisarBaseHDF.getImageData, plus two things slim GCOVs need: a
        float16 band is returned as float32, and rtcGammaToSigmaFactor is
        read from the shared factor file (<factorFrom>/<granule name>, the
        symlink downloadNISARoptimized leaves) when the granule lacks it.
        '''
        grid = self.h5[self.product][self.bands][self.frequency]
        if productField == 'rtcGammaToSigmaFactor' and productField not in grid:
            if self.factorFrom is None:
                raise ValueError(
                    'rtcGammaToSigmaFactor is not in this (slim) GCOV; pass '
                    'factorFrom=<downloadNISARoptimized --factorDir> to openHDF')
            setattr(self, productField, self._sharedFactor(useNumpy))
        else:
            super().getImageData(productField, layer=layer, useNumpy=useNumpy,
                                 power=power)
        data = getattr(self, productField, None)
        if useNumpy and isinstance(data, np.ndarray) and data.dtype == np.float16:
            setattr(self, productField, data.astype(np.float32))

    def _sharedFactor(self, useNumpy):
        '''Read the shared rtcGammaToSigmaFactor for this granule.'''
        path = os.path.join(self.factorFrom, os.path.basename(str(self.hdfFile)))
        h5path = f'/science/{self.sar}/{self.product}/{self.bands}/' \
            f'{self.frequency}/rtcGammaToSigmaFactor'
        with h5py.File(path, 'r') as f:
            data = f[h5path][()]
        gridShape = self.h5[self.product][self.bands][self.frequency][
            'mask'].shape
        # A factor on another grid (partial frame one cycle, full the next)
        # would silently misregister the gamma->sigma conversion
        if data.shape != gridShape:
            raise ValueError(f'shared factor {path} is {data.shape}, granule '
                             f'grid is {gridShape}')
        return self.downsample(data)

    def computeSigma(self):
        '''
        Convert data from gamma to sigma
        Returns
        -------
        None.

        '''
        rtcConversion = getattr(self, 'rtcGammaToSigmaFactor')
        if rtcConversion is None:
            print('Warning: No rtcGammaToSigmaFactor to do conversion')
        # check it has been converted already
        if self.backscatterType == 'sigma0':
            return 
        for field in self.dataFields:
            if field in ['HHHH', 'VVVV', 'HHVV', 'VVHH']:
                setattr(self, field, self._computeSigma(getattr(self, field),
                                                        rtcConversion))
        self.backscatterType = 'sigma0'

    def _computeSigma(self, data, rtcConversion):
        return rtcConversion * data

    def computedB(self):
        '''
        Convert data from gamma to sigma
        Returns
        -------
        None.

        '''
        if self.units == 'dB':
            return
        for field in self.dataFields:
            if field in ['HHHH', 'VVVV', 'HHVV', 'VVHH']:
                setattr(self, field, self._computedB(getattr(self, field)))
        self.units = 'dB'
                
    def _computedB(self, data):
        return 10.0 * np.log10(data)

 