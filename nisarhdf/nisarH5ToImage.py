#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 15:38:36 2025

@author: ian
"""
import argparse
import nisarhdf
from datetime import datetime
import sys
import psutil

productTypes = ['RSLC', 'ROFF', 'RIFG', 'RUNW', 'GCOV', 'GUNW', 'GOFF' ]

# Dictionary with information on other fields
fieldsDict = {'ROFF': {'pixelOffsets': ['slantRangeOffset',
                                        'slantRangeOffsetVariance',
                                        'alongTrackOffset', 
                                        'alongTrackOffsetVariance',
                                        'crossOffsetVariance', 
                                        'correlationSurfacePeak',
                                        'snr', 
                                        'digitalElevationModel']},
              'RIFG': {'interferogram': ['wrappedInterferogram',
                                         'coherenceMagnitude']},
              'RUNW': {'interferogram': ['unwrappedPhase',
                                         'coherenceMagnitude',
                                         'connectedComponents', 
                                         'ionospherePhaseScreen',
                                         'ionospherePhaseScreenUncertainty',
                                         'digitalElevationModel']},
              'RSLC': {None: ['HH', 
                              'VV',
                              'HV', 
                              'VH']},
              'GCOV': {None: ['HHHH', 
                              'VVVV',
                              'HVHV', 
                              'VHVH',
                              'mask',
                              'numberOfLooks',
                              'rtcGammaToSigmaFactor']},
              'GUNW': {'unwrappedInterferogram': ['unwrappedPhase',
                                                  'coherenceMagnitude',
                                                  'connectedComponents',
                                                  'ionospherePhaseScreen',
                                                  'ionospherePhaseScreenUncertainty'], 
                       'wrappedInterferogram': ['wrappedInterferogram',
                                                'coherenceMagnitude']}, 
              'GOFF': {'pixelOffsets': ['slantRangeOffset',
                                        'slantRangeOffsetVariance',
                                        'alongTrackOffset',
                                        'alongTrackOffsetVariance',
                                        'crossOffsetVariance', 
                                        'correlationSurfacePeak',
                                        'snr']}
              }

productTypeDict = {'ROFF': ['pixelOffsets'],
                    'RIFG': ['interferogram'],
                    'RUNW': ['interferogram'],
                    'RSLC': [None],
                    'GCOV':[None],
                    'GUNW': ['unwrappedInterferogram', 'wrappedInterferogram'],  
                    'GOFF': ['pixelOffsets']
                    }

defaultFieldsDict = {'RSLC': {None: ['HH', 'VV', 'HV', 'VH']},
                     'ROFF': {'pixelOffsets': ['slantRangeOffset', 
                                               'alongTrackOffset']},
                     'RIFG': {'interferogram': ['wrappedInterferogram']},
                     'RUNW': {'interferogram': ['unwrappedPhase']},
                     'GCOV': {None: ['HHHH', 'VVVV', 'HVHV', 'VHVH']},
                     'GUNW': {'unwrappedInterferogram': ['unwrappedPhase'], 
                              'wrappedInterferogram': ['wrappedInterferogram']}, 
                     'GOFF': {'pixelOffsets': ['slantRangeOffset',
                                               'alongTrackOffset']}
                     
                    }

quickLookDefaultDownsampleFactors = {'RSLC': 1,
                                     'ROFF': 1, # not used
                                     'RIFG': 8,
                                     'RUNW': 2,
                                     'GCOV': 8,
                                     'GUNW': 2, 
                                     'GOFF': 2} # not use
    
def specialKeywords(myArgs, applicableProducts, specialArgs, args, 
                    default=False):
    '''
    Setup product specific arguments

    Parameters
    ----------
    productType : TYPE
        DESCRIPTION.
    applicableProducts : TYPE
        DESCRIPTION.
    keyword : TYPE
        DESCRIPTION.
    args : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    for arg in specialArgs:
        if myArgs['product'] in applicableProducts:
            myArgs[arg] = getattr(args, arg)
        else:
            myArgs[arg] = default
            #print(f'--{arg} only for {applicableProducts}, ignoring')
            

def myerror(message,myLogger=None):
    """ print error and exit """
    print('\n\t\033[1;31m *** ',message,' *** \033[0m\n') 
    if myLogger != None :
        myLogger.logError(message)
    sys.exit()

def parseCommandLine():
    '''
    Handle command line args
    '''
    parser = argparse.ArgumentParser(
        description='\n\n\033[1mConvert NISAR H5 product to Tiff(s), COG(s) or binary files'
        '\033[0m\n\n',)
    parser.add_argument('productName', type=str,
                        help='NISAR h5 product (file path, s3, or https)')
    parser.add_argument('output', type=str, nargs="?",
                        help='Root name for ouput. Band (e.g., .unwrappedPhase) and file type (e.g., .tiff) will be appended. (Not required for --info)')
    parser.add_argument('--info', action="store_true",
                        help='Print summary info only')
    parser.add_argument('--quickLook', action="store_true",
                        help='Write a quick look PNG file (GCOV )')
    #parser.add_argument('--ros3', action="store_true",
    #                    help='Use ros3 for network (s3,https), slower, but more memory efficient')  
    parser.add_argument('--productFamily', type=str, default='',
                        help='NISAR product type (e.g., RUNW, GUNW etc) '
                        '[parse from product path]',
                        choices=productTypes)
    parser.add_argument('--frequencyB', action="store_true",
                        help='Select frequencyB [frequencyA] ')    
    parser.add_argument('--dB', action="store_true",
                        help='Output results in dB (GCOV only)')  
    parser.add_argument('--sigma0', action="store_true",
                        help='Output results as sigma0 (GCOV only) [gamma0]')  
    parser.add_argument('--scaleToPixels', action="store_true",
                        help='Scale offsets to pixels (ROFF, GOFF only)')  
    parser.add_argument('--wrapped', action="store_true",
                        help='Wrapped interferogram (GUNW only)')  
    parser.add_argument('--polarization', type=str, default=None,
                        choices=['HH', 'VV', 'HV', 'VH'],
                        help='Polarization for \n\n\033[1mnon-GCOV\n\n\033[0m '
                        'products [first like pol]')
    parser.add_argument('--downsampleFactor', type=int, default=None,
                        choices=[2, 4, 8],
                        help='Reduce resolution by 1/subsampleFactor (GCOV only) [full res]')
    parser.add_argument('--fields', type=str, default=[], nargs='+',                       
                        help='Select fields including GCOV covariance terms (e.g., '
                        'HHHH, HVHV), use all for everything. '
                        'Use --info to see defaults for product type [defaults]')
    parser.add_argument('--layers', type=str, default=['1', '2', '3'], nargs='+',
                        choices=['1', '2', '3'],                    
                        help='For offsets only, layer numbers to access [1-3]')
    parser.add_argument('--outputFormat', type=str, default='COG',
                        choices=['COG', 'GTiff', 'binary'],
                        help='Output format [COG: cloud optimized geotiff]')
    #
    args = parser.parse_args()
    #
    myArgs = {}
    #
    if args.output is None and not args.info:
        myerror('Unless --info set, specify an output file name')
    #
    myArgs['product'] = getProductType(args.productFamily, 
                                       args.productName)
    #
    specialKeywords(myArgs, ['GCOV'], ['dB', 'sigma0'], args)
    specialKeywords(myArgs, ['ROFF', 'GOFF'], ['scaleToPixels'], args)
    specialKeywords(myArgs, ['GUNW'], ['wrapped'], args)
    myArgs['frequency'] = {False: 'frequencyA',
                           True: 'frequencyB'}[args.frequencyB]
    #
    processHDFOpenKeywords(args, myArgs)
    #
    myArgs['productType'] = productTypeDict[myArgs['product']][0]
    if myArgs['wrapped']:
        myArgs['productType'] = productTypeDict[myArgs['product']][1]
    #
    myArgs['fields'] = reconcileFields(args, myArgs)
    #
    # args that need no checking
    for arg in ['productName', 'output', 'polarization',
                'outputFormat', 'info',
                'downsampleFactor', 'quickLook']:
        myArgs[arg] = getattr(args, arg)
    myArgs['ros3'] = True
    #
    if myArgs['quickLook'] and myArgs['product'] in ['ROFF', 'GOFF']:
        print(f'quickLooks not supported for ROFF and GOFF')
        myArgs['quickLook'] = False
    elif myArgs['quickLook']:
        if myArgs['downsampleFactor'] == 0:
             myArgs['downsampleFactor'] = quickLookDefaultDownsampleFactors[myArgs['product']]
             if myArgs['product'] in ['GUNW'] and myArgs['wrapped']:
                 myArgs['downsampleFactor'] *= 4
    if myArgs['downsampleFactor'] is None:
        # If not set by one of the defaults, default to 1
        myArgs['downsampleFactor'] = 1
        
    #
    return myArgs

def  processHDFOpenKeywords(args, myArgs):
    '''
    Process the set of keywords that are passed to the hdfOpen method
    '''
    myArgs['hdfOpenKeywords'] = {}
    # keywords that only apply to GCOV
    if myArgs['product'] in ['GCOV']:
        myArgs['hdfOpenKeywords']['dB'] = myArgs['dB']
        myArgs['hdfOpenKeywords']['sigma0'] = myArgs['sigma0']
        if args.polarization is not None:
            myerror('--polarization is not valid for GCOV. \n'
                    'Polarization specified through fields (e.g --fields HHHH HVHV)')
    else:
        if myArgs['product'] not in ['RSLC']:
            myArgs['hdfOpenKeywords']['polarization'] = args.polarization
    #
    # Only offset products have layers
    if myArgs['product'] in ['GOFF', 'ROFF']:
        myArgs['hdfOpenKeywords']['layers'] = [f'layer{x}' for x in args.layers]


def reconcileFields(args, myArgs):
    '''
    Use defaults and command line options to produce final set of fields
    '''
    #
    # Reconcile fields
    fields = args.fields
    product = myArgs['product']
    if len(fields) == 0:
        fields = defaultFieldsDict[product][myArgs['productType'] ]
    elif 'all' in fields:
        fields = fieldsDict[product][myArgs['productType']]
    else:
        for field in fields:
            if field not in fieldsDict[product][myArgs['productType']]:
                myerror(f'{field} not in '
                        f'{fieldsDict[product][myArgs["productType"]] }') 
    return fields

def printError(msg):
        '''
        Print error message
        Parameters
        ----------
        msg : str
            error message.
        Returns
        -------
        None

        '''
        print(f'\n\t\033[1;31m *** {msg} *** \033[0m\n')


def getProductType(productType, productPath):
    '''
    Return either passed in productType if or find type from productPath
    '''
    if productType in productTypes:
        return productType
    # Find product Type
    for productType in productTypes:
        if productType in productPath.upper():
            return productType
    printError('Product type not specified and cannot parse from product path')
    
def outputData(myArgs, myProduct):
    '''
    Output data

    Parameters
    ----------
    myArgs : dict
        Parameters.
    myProduct : nisarhdf obj
        The product.

    Returns
    -------
    None.

    '''
    keywords = {}
    if myArgs['productType'] in ['GOFF', 'ROFF']:
        keywords = {'scaleToPixels': myArgs['scaleToPixels']}
    #if myArgs['product'] in ['GCOV']: 
     #   keywords = {'dB': myArgs['dB'], 'sigma0': myArgs['sigma0']}
        #print(keywords)
    #
    #for myArgs[]
    tiff, driverName = {'GTiff': [True, 'GTiff'],
                        'COG': [True, 'COG'],
                        'binary': [False, '']}[myArgs['outputFormat']]
    # write data 
    myProduct.writeData(myArgs['output'],
                        tiff=tiff,
                        quickLook=myArgs['quickLook'],           
                        driverName=driverName, bands=None, **keywords)
    return
 

def printInfo(myArgs, myProduct):
    '''

    Parameters
    ----------
    myArgs : dict
        CLI parameters.
    myProduct : nisarhdf object
        The current product.
    Returns
    -------
    None.

    '''
    #
    # List defaults
    if myArgs["product"] in ['GCOV', 'RSLC']:
        print("\n\033[1mDefault polarization:\033[0m")
        for pol in myProduct.polarizations:
            print(pol, end=' ')
        print('')
        print('\033[1mAvailable fields: \033[0m', end='\n')
        if myArgs["product"] in ['GCOV']:
            prodSpecTerms = myProduct.covTerms
        else:
            prodSpecTerms = myProduct.polarizations
        for pol in prodSpecTerms:
            print(f'{pol}')
        for field in fieldsDict[myArgs["product"]][None]:
            if field[0:2] not in ['HH', 'VV', 'HV', 'VH']:
                print(f'{field}')
    else: 
        for productType in defaultFieldsDict[myArgs["product"]]:
            prodName = {True: myArgs["product"],
                        False: myArgs["productType"]}[productType is None]
            print(f'\033[1m\nDefaults for {prodName}: \033[0m', end='\n')
            for field in defaultFieldsDict[myArgs["product"]][productType]:
                print(f'{field}')
            print(f'\033[1mAvailable fields for {productType}: \033[0m', end='\n')
            for field in fieldsDict[myArgs["product"]][productType]:
                print(f'{field}')
    #
    myProduct.printParams()


def run():
    ''' 
    Main driver for program
    '''
    myArgs = parseCommandLine()
    net_before = psutil.net_io_counters()
    # Debug
    if False:
        for key in myArgs:
            print(key, myArgs[key])
    #sys.exit()
    #
    #
    myProduct = getattr(nisarhdf,
                        f'nisar{myArgs["product"]}HDF')(
                            frequency=myArgs['frequency'])
    # read data
    start = datetime.now()
    myProduct.openHDF(myArgs['productName'],
                      productType=myArgs['productType'],
                      noLoadData=myArgs['info'],
                      fields=myArgs['fields'],
                      useRos3=myArgs['ros3'],
                      page_buf_size=1 * 1024**3,
                      downsampleFactor=myArgs['downsampleFactor'],
                      **myArgs['hdfOpenKeywords']
                      )
    read = datetime.now()
    #
    # print info and return for --info
    if myArgs['info']:
        printInfo(myArgs, myProduct)
    else:
        outputData(myArgs, myProduct)  
    write = datetime.now()
    print('\n\033[1mRun time: \033[0m', end='\n')
    print('Load data', read-start)
    print('Write data', write-read)
    print('Total', write-start)
    #
    print('\n\033[1mNetwork traffic: \033[0m', end='\n')
    net_after = psutil.net_io_counters()    
    bytes_sent = net_after.bytes_sent - net_before.bytes_sent
    bytes_recv = net_after.bytes_recv - net_before.bytes_recv
    print(f"Bytes sent: {bytes_sent/1e6:.3f} MB")
    print(f"Bytes received: {bytes_recv/1e6:.3f} MB")

if __name__ == "__main__":
    run()
    
