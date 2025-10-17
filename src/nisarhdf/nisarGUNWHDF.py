#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 15:21:35 2024

@author: ian
"""

from nisarhdf import nisarBaseGeocodedHDF


class nisarGUNWHDF(nisarBaseGeocodedHDF):
    '''
    This class creates objects to work with nisar RUNWHDF images.
    '''

    def __init__(self,  sar='LSAR', product='GUNW', frequency='frequencyA',
                 productType='unwrappedInterferogram', polarization=None,
                 layer=None, productData='unwrappedPhase', bands='grids',
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
        self.productType = productType
        self.productParams = ['NumberRangeLooks', 'NumberAzimuthLooks']
        self.lookType = None

    def parseParams(self, fields=None, secondary=False, noLoadData=False,
                    productType='unwrappedInterferogram', polarization=None,
                    **keywords):
        '''
        Parse all the params needed to make a geodatNRxNA.geojson file

        Returns
        -------
        None.

        '''
        self.productType = productType
        if not secondary:
            self.getPolarization(polarization)
        #
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
        self.getInterferogramPixelOffsets()
        #
        fieldDict = {'unwrappedInterferogram':
                     ['coherenceMagnitude',
                      'connectedComponents',
                      'ionospherePhaseScreen',
                      'ionospherePhaseScreenUncertainty',
                      'unwrappedPhase'],
                     'wrappedInterferogram':
                     ['coherenceMagnitude',
                      'wrappedInterferogram']}
        if fields is None:
            fields = fieldDict[self.productType]
        self.loadData(fields, noLoadData=noLoadData)
