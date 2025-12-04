#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 08:05:24 2024

@author: ianf
"""

from nisarhdf import nisarBaseRangeDopplerHDF
import os
import numpy as np


class nisarRSLCHDF(nisarBaseRangeDopplerHDF):
    '''
    This class creates objects to work with nisar RUNWHDF images.
    '''

    def __init__(self,  sar='LSAR', frequency='frequencyA', product='RSLC',
                 polarization='HH', byteOrder='LSB', productType=None,
                 bands='swaths', productData=None):
        '''
        NISAR RSLC HDF Reader

        Parameters
        ----------
        sar : str, optional
            SAR Idendtifier (always LSAR for now). The default is 'LSAR'.
        frequency : str, optional
            frequency band to extract. The default is 'frequencyA'.
        polarization : str, optional
            Polarization. The default is 'HH'.

        Returns
        -------
        None.

        '''
        nisarBaseRangeDopplerHDF.__init__(self,
                                          sar=sar,
                                          product=product,
                                          frequency=frequency,
                                          productType=productType,
                                          polarization=polarization,
                                          # layer=layer,
                                          productData=productData,
                                          bands=bands,
                                          bytOrder=byteOrder)
        self.lookType = 'SLC'
        self.productParams = ['NumberRangeLooks', 'NumberAzimuthLooks']
        for param in self.RDParams:
            self.productParams.append(f'{self.lookType}{param}')
        self.isSecondary = False

    def parseParams(self,
                    noLoadData=False,
                    fields=None,
                    downsampleFactor={'downsampleFactorRow': 1,
                                      'downsampleFactorColumn': 1},
                    useNumpy=False,
                    power=False,
                    **keywords):
        '''
        Parse all the params needed for offsets

        Returns
        -------
        None.

        '''

        # self.getOffsetParams()
        # self.getNumberOfLooks()
        if 'polarization' not in keywords:
            polarization = None
        self.getPolarization(polarization)
        self.getOrbitAndFrame(**keywords)
        self.getLookDirection()
        self.getOrbitPassDirection()
        self.parseRefDate()
        self.getSLCSLantRangeAndZeroDoppler()
        self.effectivePRF()
        self.orbit = \
            self.parseStateVectors(XMLOrbit=self.referenceOrbitXML)               
        self.getSceneCenterSatelliteHeight()
        self.getSLCSize(**keywords)
        self.NumberRangeLooks = 1
        self.NumberAzimuthLooks = 1
        #
        self.getPolarizations()
        if fields is None:    
            fields = self.polarizations
        fields = [x for x in fields if x in self.polarizations]
        # load data
        print('loading as power', power)
        
        self.loadData(fields, useNumpy=useNumpy, power=power,
                      noLoadData=noLoadData)
        self.NumberRangeLooks *= self.downsampleFactorColumn
        self.NumberAzimuthLooks *= self.downsampleFactorRow
        #  
        self.getRangeBandWidth()
        self.getMLSize()
        self.getMLZeroDopplerTime(SLC=True)
        self.getCorrectedTime()
        self.getMLSlantRange()
        self.getCorners()
        self.getCenterLatLon()
        self.getRangeErrorCorrection()
        self.getTimeToFirstSample()
        self.getSkewOffsets()
        #myRSLC.getCenterIncidenceAngle()
        self.getSquint()
        self.getDeltaT()
        self.getCenterLatLon()
        
        self.getSceneCenterSatelliteHeight()
        self.getWavelength()
        _, self.MLIncidenceCenter = self.computeAngles(self.MLCenterRange,
                                                            self.MLMidZeroDopplerTime, 
                                                            np.array([0]),
                                                            degrees=False)
        _, self.SLCIncidenceCenter = self.computeAngles(self.SLCCenterRange,
                                                            self.SLCMidZeroDopplerTime, 
                                                            np.array([0]),
                                                            degrees=False)
        
        #
        # if the product has been multilooked, add ML params
        if self.NumberRangeLooks > 1 or self.NumberAzimuthLooks > 1:
            for param in self.RDParams:
                self.productParams.append(f'ML{param}')


    
        
    def writeSLC(self, fileName, frequency='frequencyA', dataType='>c8'):
        '''
        Parameters
        ----------
        fileName : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        swaths = self.h5Full['science'][self.sar][self.product]['swaths']
        data = np.array(swaths['frequencyA']['HH'].astype('c8'))
        self._writeImageData(fileName, data, dataType, noData=-2.e9)

    def writeOldParFile(self, filename):
        '''
        Write an old gamma style par file for SLC. GrIMP only.
        '''
        with open(filename, 'w') as fp:
            print(os.path.basename(self.hdfFile), file=fp)
            print(self.datetime.strftime('%d %b %Y').upper(), file=fp)
            print(self.datetime.strftime('%H %M %S.%f'), file=fp)
            print(' ', file=fp)
            print('0. 0.\n0.', file=fp)
            print(self.SpaceCraftAltitude, file=fp)
            print('0.0\n 0.0 0.0 0.0\n 0.0 0.0 0.0\n 0.0 0.0 0.0', file=fp)
            print(self.PRF, file=fp)
            print('0.0 0.0 0.0 0.0\n0.\n0.\n0. 0. 0.', file=fp)
            print(f'{self.SLCNearRange:0.2f}  {self.SLCCenterRange:.2f} '
                  f'{self.SLCFarRange:.2f}', file=fp)
            print(f'{self.SLCRangePixelSize:.8f}  0', file=fp)
            print('0\n1\n0\n0\n0\n1\n1 1\n0.', file=fp)
            print(f'{self.SLCAzimuthPixelSize:.8f} 0.', file=fp)
            print(f'{self.SLCRangeSize} {self.SLCAzimuthSize}', file=fp)
            print(f'{self.orbit.StateVectorInterval}', file=fp)
            print(f'{self.orbit.TimeOfFirstStateVector}', file=fp)
            for p, v in zip(self.orbit.position, self.orbit.velocity):
                print(f'{p[0]:10e} {p[1]:10e} {p[2]:10e}', file=fp)
                print(f'{v[0]:e} {v[1]:e} {v[2]:e}', file=fp)


    def getSLCSize(self, frequency='frequencyA', **kwds):
        '''
        Get size of SLC image

        Returns
        -------
        None.

        '''
        swath = self.h5[self.product]['swaths']
        self.SLCAzimuthSize, self.SLCRangeSize = \
            swath[frequency][self.polarization].shape

    def getDateTime(self):
        '''
        Get the time
        Returns
        -------
        None.

        '''
        id = self.h5Full['science'][self.sar]['identification']
        dateStr = np.array(id["zeroDopplerStartTime"]).item().decode('ascii')
        self.datetime, self.Date = self.parseDateStr(dateStr)

    def getSLCSLantRangeAndZeroDoppler(self, frequency='frequencyA'):
        '''
        Get near, center, and far slant range

        Returns
        -------
        None.

        '''
        swaths = self.h5[self.product]['swaths']
        self.SLCNearRange = swaths[frequency]['slantRange'][0]
        self.SLCFarRange = swaths[frequency]['slantRange'][-1]
        self.SLCCenterRange = (self.SLCNearRange + self.SLCFarRange) * 0.5
        self.SLCRangePixelSize = swaths[frequency]['slantRange'][1] - \
            swaths[frequency]['slantRange'][0]
        self.SLCRangePixelSize = np.array(
            swaths[frequency]['slantRangeSpacing']).item()
        #
        self.SLCFirstZeroDopplerTime = swaths['zeroDopplerTime'][0]
        self.SLCLastZeroDopplerTime = swaths['zeroDopplerTime'][-1]
        self.SLCMidZeroDopplerTime = \
            (self.SLCFirstZeroDopplerTime + self.SLCLastZeroDopplerTime) * 0.5
        self.SLCZeroDopplerTimeSpacing = np.array(swaths['zeroDopplerTimeSpacing']).item()
        self.PRF = 1. / self.SLCZeroDopplerTimeSpacing
        #
        self.SLCSceneCenterAlongTrackSpacing()