#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 10:44:13 2024

@author: ian
"""
from nisarhdf.nisarBaseRangeDopplerHDF import nisarBaseRangeDopplerHDF
import os
import numpy as np

from nisarhdf import writeMultiBandVrt, formatGeojson

class nisarRUNWHDF(nisarBaseRangeDopplerHDF):
    '''
    This class creates objects to work with nisar RUNWHDF images.
    '''

    def __init__(self,  sar='LSAR', frequency='frequencyA',
                 polarization=None, isSecondary=False, 
                 productType='interferogram',
                 referenceOrbitXML=None, secondaryOrbitXML=None, debug=False):
        '''
       sar : str, optional
            SAR Idendtifier (always LSAR for now). The default is 'LSAR'.
        frequency : str, optional
            frequency band to extract. The default is 'frequencyA'.
        polarization : str, optional
            Polarization. The default is 'HH'.
        isSecondary : Boolean, optional
            For internal use only. The default is False.
        referenceOrbitXML : str, optional
            XML file to override orbit in hdf The default is None.
        secondaryOrbitXML : TYPE, optional
            XML file for secondary orbit. The default is None

        Returns
        -------
        None.

        '''
        print('ref orbit', referenceOrbitXML)
        nisarBaseRangeDopplerHDF.__init__(self,
                                          sar=sar,
                                          product='RUNW',
                                          frequency=frequency,
                                          productType=productType,
                                          polarization=polarization,
                                          layer=None,
                                          productData='unwrappedPhase',
                                          bands='swaths',
                                          isSecondary=isSecondary,
                                          referenceOrbitXML=referenceOrbitXML,
                                          secondaryOrbitXML=secondaryOrbitXML,
                                          debug=debug)
        self.lookType = 'ML'
        self.productParams = ['NumberRangeLooks', 'NumberAzimuthLooks']
        for param in self.RDParams:
            self.productParams.append(f'{self.lookType}{param}')

    def parseParams(self, fields=None, productType='interferogram', 
                    polarization=None, secondary=False, 
                    noLoadData=False, **keywords):
        '''
        Parse all the params needed to make a geodatNRxNA.geojson file

        Returns
        -------
        None.

        '''
        if not secondary:
            self.getPolarization(polarization)
        self.getGranuleNames()
        self.getOrbitAndFrame(**keywords)
        self.getLookDirection()
        self.parseRefDate()
        self.getNumberOfLooks(prodType='interferogram')
        self.getSLCSlantRange()
        self.getSLCZeroDopplerTime(secondary=secondary)
        self.getMLSize()
        self.getMLSlantRange()
        self.effectivePRF()
        self.getRangeErrorCorrection()
        if not secondary:
            self.getInterferogramPixelOffsets()
            self.getMLZeroDopplerTime()
        else:
            self.ImageName = 'secondary'
            self.getMLZeroDopplerTime(secondary=True)
            #
        self.getCorrectedTime()
        # Note if secondary, it will have been passed the reference (see below)
        self.orbit = self.parseStateVectors(XMLOrbit=self.referenceOrbitXML)
        #
        self.getOrbitPassDirection()
        self.getCorners()
        self.getRangeErrorCorrection()
        self.getTimeToFirstSample()
        self.getSkewOffsets()
        self.getCenterIncidenceAngle()
        self.getSquint()
        self.getDeltaT()
        self.getCenterLatLon()
        self.getSceneCenterSatelliteHeight()
        #self.getOrbitPassDirection()
        self.getWavelength()
        #self.getSingleLookPixelSize()
        self.getExtent()
        #
        # If not secondary, create the secondary and parse
        if not secondary:
            self.secondary = \
                nisarRUNWHDF(isSecondary=True,
                             referenceOrbitXML=self.secondaryOrbitXML,
                             debug=self.debug)
            self.secondary.h5 = self.h5
            self.secondary.parseParams(secondary=True, referenceOrbit=self.secondaryOrbit)
            #self.secondaryDatetime = self.secondary.datetime
            #self.secondaryDate = self.secondary.datetime
            self.dT = np.round((self.secondaryDatetime -
                               self.datetime).total_seconds()/86400.,
                               decimals=3)
            self.secondary.referenceOrbit = self.secondaryOrbit
            self.secondary.frame = self.frame
            self.genGeodatProperties()
            if fields is None:
                fields = ['coherenceMagnitude',
                          'connectedComponents',
                          'ionospherePhaseScreen',
                          'ionospherePhaseScreenUncertainty',
                          'unwrappedPhase',
                          'digitalElevationModel']
            self.loadData(fields, noLoadData=noLoadData)

