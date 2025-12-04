import argparse
import nisarhdf
from datetime import datetime
import sys
import psutil
import os

def positive_int(value):
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError(f"{value} must be an integer > 0")
    return ivalue

def parseCommandLine():
    '''
    Handle command line args
    '''
    parser = argparse.ArgumentParser(
        description='\n\n\033[1mMultilook NISAR RSLC to power image Tiff(s), COG(s) or binary files optionally to dB'
        '\033[0m\n\n',)
    parser.add_argument('RSLCFile', type=str,
                        help='NISAR RSLC product (file path, s3, or https)')
    parser.add_argument('output', type=str, nargs="?",
                        help='Root name for ouput. file type (e.g., .tiff) will be appended.')
    parser.add_argument('--frequencyB', action="store_true",
                        help='Select frequencyB [frequencyA] ')    
    parser.add_argument('--dB', action="store_true",
                        help='Output results in dB (GCOV only)')  
    parser.add_argument('--polarization', type=str, default='HH',
                        choices=['HH', 'VV', 'HV', 'VH'],
                        help='Polarization for \n\n\033[1mnon-GCOV\n\n\033[0m '
                        'products [first like pol]')
    parser.add_argument('--orbitPOE', type=str, default=None,
                        help='Precision Orbit file [orbit from RSLC]')
    parser.add_argument('--nlooksRange', type=positive_int, default=1,
                        help='Number of looks in range [1]')
    parser.add_argument('--nlooksAzimuth', type=positive_int, default=1,
                        help='Number of looks in azimuth [1]')
    parser.add_argument('--outputFormat', type=str, default='COG',
                        choices=['COG', 'GTiff', 'binary'],
                        help='Output format [GTiff: geotiff]')
    parser.add_argument('--geojsonFile', type=str, default=None,
                        help='Override default name [geodatnlooksRangexnlooksAzimuth.geojson]')
    parser.add_argument('--noGeojson', action="store_true",
                        help='Do not output a geojson file')  
    parser.add_argument('--geojsonOnly', action="store_true",
                        help='Only create a geojson file with no image')  
    #
    args = parser.parse_args()
    #
    myArgs = {}
    #
    for arg in ['RSLCFile', 'output', 'polarization',
                'outputFormat', 'dB', 'orbitPOE', 'geojsonOnly']:
        myArgs[arg] = getattr(args, arg)
    myArgs['outputPath'] = os.path.dirname(myArgs['output'])
    #if args.output is None and not args.info:
    #    myerror('Unless --info set, specify an output file name')
    #
    myArgs['frequency'] = {False: 'frequencyA',
                           True: 'frequencyB'}[args.frequencyB]
    #
    if args.geojsonFile is None and not args.noGeojson:
        geojsonFile  = f'geodat{args.nlooksRange}x{args.nlooksAzimuth}.geojson'
        myArgs['geojsonFile'] = os.path.join(myArgs['outputPath'], geojsonFile)
    else:
        myArgs['geojsonFile'] = args.geojsonFile
    #
    myArgs['downsampleFactor'] = {'downsampleFactorRow': args.nlooksAzimuth,
                                  'downsampleFactorColumn': args.nlooksRange}
    #
    return myArgs


def print_s3_files_flat(tree, nameFilter=None, fileEndsWith=None, excludeName=None):
    """
    Print full paths of all files in the nested S3 dictionary (flat list),
    with optional filters.

    Parameters
    ----------
    tree : dict
        Nested dict returned by list_s3_tree.
    nameFilter : str or None
        Include only names containing this substring (applies to dirs and files)
    fileEndsWith : str or None
        Include only files ending with this string
    excludeName : str or None
        Skip names containing this substring (applies to dirs and files)
    """
    # Print files at this level
    for f in tree.get("files", []):
        fname = f.split("/")[-1]
        if excludeName and excludeName in fname:
            continue
        if nameFilter and nameFilter not in fname:
            continue
        if fileEndsWith and not fname.endswith(fileEndsWith):
            continue
        print(f)

    # Recurse into subdirectories
    for dirname, subtree in tree.get("directories", {}).items():
        if excludeName and excludeName in dirname:
            continue
        if nameFilter and nameFilter not in dirname:
            continue
        print_s3_files_flat(subtree, nameFilter=nameFilter,
                            fileEndsWith=fileEndsWith, excludeName=excludeName)

def run():
    ''' 
    Main driver for program
    '''
    myArgs = parseCommandLine()
    start = datetime.now()
    net_before = psutil.net_io_counters()
    # Debug
    if True:
        for key in myArgs:
            print(key, myArgs[key])
    #sys.exit()
    # Read SLC
    myRSLC = nisarhdf.nisarRSLCHDF()
    noLoadData = myArgs['geojsonOnly'] # Use for debug
    myRSLC.openHDF(myArgs['RSLCFile'],
                 fields=[myArgs['polarization']], 
                 useRos3=False,
                 useNumpy=True,
                 power=True,
                 downsampleFactor=myArgs['downsampleFactor'],
                 referenceOrbitXML=myArgs['orbitPOE'],
                 noLoadData=noLoadData)
    #
    read = datetime.now()
    #
    if myArgs['outputFormat'] == 'binary':
        tiff = False
        driverName = None
    elif myArgs['outputFormat'] in ['GTiff', 'COG']:
        driverName = myArgs['outputFormat']
        tiff = True
    #
    if not noLoadData:
        myRSLC.writeData(myArgs['output'],
                         bands=[myArgs['polarization']],
                         tiff=tiff,
                         grimp=True,
                         driverName=driverName,
                         vrtFile=f"{myArgs['output']}.vrt",
                         noSuffix=True
                         )
    #
    if myArgs['geojsonFile'] is not None:
        myRSLC.writeGeodatGeojson(filename=myArgs['geojsonFile'])
    #
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
    