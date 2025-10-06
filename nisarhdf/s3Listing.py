
import argparse
import subprocess
import sys
from datetime import datetime
from urllib.parse import urlparse

s3Count = [0]

def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='\n\n\033[1mGet S3 recursive bucket listing\033[0m\n\n',)
    parser.add_argument('s3link', type=str,
                        help='s3 bucket to query')
    parser.add_argument('--nameFilter', type=str, default=None,
                    help='Only include paths containing this substring')
    parser.add_argument('--excludeName', type=str, default=None,
                    help='Exclude paths containing this substring')
    parser.add_argument('--fileEndsWith', type=str, default=None,
                    help='Only include files ending with this string, e.g., ".h5"')
    parser.add_argument('--maxLevels', type=int, default=None,
                    help='Maximum recursion depth for directories; None = recurse fully')
    parser.add_argument('--createdAfter', type=str, default=None,
                    help='Only show files created after specified data YYYY-MM-DD')
    # Flat flag
    parser.add_argument('--flat', action='store_true',
                    help='Print files as a flat list of full paths instead of a hierarchy')
    # Long flag
    parser.add_argument('--long', action='store_true',
                    help='Print all file information (date, size)')
    args = parser.parse_args()
    printKeywords, searchKeywords = {}, {}
    for arg in {'nameFilter', 'fileEndsWith', 'excludeName', 'long'}:
        printKeywords[arg] = getattr(args, arg)
    if args.createdAfter is not None:
        printKeywords['createdAfter'] = datetime.strptime(args.createdAfter,
                                                          "%Y-%m-%d")
    else:
        printKeywords['createdAfter'] = None
    #
    #for arg in {'maxLevels'}:
    #    searchKeywords[arg] = getattr(args, arg)
    #
    return args.s3link, searchKeywords, printKeywords, args.flat
    
def human_size(num_bytes):
    '''
    Convert num_bytes to format like ls -h
    '''
    for unit in ["B","K","M","G","T","P","E"]:
        if abs(num_bytes) < 1024.0:
            return f"{num_bytes:6.1f}{unit}"
        num_bytes /= 1024.0
    return f"{num_bytes:7.1f}Z"  # fall back for huge numbers

def is_s3_dir(s3_uri, profile=None):
    '''
    Make sure is directory. Add trailing slash if it is.
    '''
    # Ensure command
    cmd = ["aws", "s3", "ls"]
    if profile:
        cmd += ["--profile", profile]
    # If already ends with /, it's a directory
    if s3_uri.endswith("/"):
        return True
    # Try listing with trailing slash
    try:
        subprocess.check_output(cmd + [s3_uri + "/"], text=True)
        return True
    except subprocess.CalledProcessError:
        return False

def normalize_key(s3_uri, key):
    '''
    Assemble full path, avoiding any redudancies.
    '''
    # Parse bucket + prefix
    parsed = urlparse(s3_uri)
    bucket = parsed.netloc
    prefix = parsed.path.lstrip("/")   # remove leading /

    # Ensure prefix ends with /
    if prefix and not prefix.endswith("/"):
        prefix += "/"

    # If key starts with prefix, strip it
    if key.startswith(prefix):
        key = key[len(prefix):]

    return f"s3://{bucket}/{prefix}{key}"

    
def list_s3_tree_recursive(s3_uri, profile=None):
    """
    Build a nested dictionary of S3 contents.
    Handles both directories/prefixes and single files.
    Uses one --recursive call for directories.
    """
    tree = {"files": [], "info": [], "date": [], "directories": {}}

    cmd_base = ["aws", "s3", "ls"]
    if profile:
        cmd_base += ["--profile", profile]

    # Check if URI is a prefix/directory (ends with / or is listable)
    is_dir = is_s3_dir(s3_uri, profile=profile)
    if is_dir:
        # Ensure trailing slash
        if not s3_uri.endswith("/"):
            s3_uri += "/"

        # Do a single recursive ls
        cmd = cmd_base + ["--recursive", s3_uri]
        try:
            result = subprocess.check_output(cmd, text=True)
        except subprocess.CalledProcessError:
            print(f"\033[1;31mError: Check if {s3_uri} exists\033[0m")
            return tree

        for line in result.splitlines():
            parts = line.split()
            if len(parts) >= 4:
                date_str = parts[0]
                size_str = parts[2]
                key = " ".join(parts[3:])
                full_key = normalize_key(s3_uri, key)
                size = int(size_str)

                path_parts = key.strip("/").split("/")

                current = tree
                for folder in path_parts[:-1]:
                    folder += "/"
                    if folder not in current["directories"]:
                        current["directories"][folder] = {
                            "files": [],
                            "info": [],
                            "date": [],
                            "directories": {}
                        }
                    current = current["directories"][folder]
                current["files"].append(full_key)
                current["info"].append(f"{date_str} {human_size(size)}")
                current["date"].append(datetime.strptime(date_str, "%Y-%m-%d"))
    else:
        # Single file
        try:
            result = subprocess.check_output(cmd_base + [s3_uri], text=True)
            for line in result.splitlines():
                parts = line.split()
                if len(parts) >= 4:
                    date_str = parts[0]
                    size_str = parts[2]
                    key = " ".join(parts[3:])
                    size = int(size_str)
                    tree["files"].append(s3_uri)
                    tree["info"].append(f"{date_str} {human_size(size)}")
                    tree["date"].append(datetime.strptime(date_str, "%Y-%m-%d"))
        except subprocess.CalledProcessError:
            print(f"\033[1;31mError: File {s3_uri} not found\033[0m")

    return tree


def print_s3_tree(tree, name="", level=0, 
                  nameFilter=None, fileEndsWith=None, excludeName=None,
                  long=False, createdAfter=None):
    """
    Pretty-print the nested S3 dictionary with indentation.
    Supports filters:
      - nameFilter: include only names containing this substring (applies to dirs and files)
      - fileEndsWith: include only files ending with this string
      - excludeName: skip names containing this substring (applies to dirs and files)

    Parameters
    ----------
    tree : dict
        Nested dict returned by list_s3_tree.
    name : str
        Name of the current directory (optional).
    level : int
        Indentation level.
    long : bool
        Print all file info.
    createdAfter: datestr
        Only return files with creation dates earlier than given date. 
    """
    indent = " " * level

    # Check if this directory should be printed
    if name:
        if excludeName and excludeName in name:
            show_dir = False
        elif nameFilter is None or nameFilter in name:
            print(f"{indent}{name}/")
            show_dir = True
        else:
            show_dir = False
    else:
        show_dir = True  # top-level bucket

    # Print files
    for f, info, date in zip(tree.get("files", []),
                             tree.get("info", []),
                             tree.get("date", [])):
        fname = f.split("/")[-1]
        if createdAfter is not None and date < createdAfter:
            continue
        if excludeName and excludeName in fname:
            continue
        if nameFilter and nameFilter not in fname:
            continue
        if fileEndsWith and not fname.endswith(fileEndsWith):
            continue
        if not long:
            info = ''
        print(f"{indent} {info} {fname}")

    # Recurse into subdirectories
    for dirname, subtree in tree.get("directories", {}).items():
        #if excludeName and excludeName in dirname:
         #   continue
        #if nameFilter and nameFilter not in dirname:
        #    continue
        print_s3_tree(subtree, dirname.rstrip("/"), level+1,
                      nameFilter=nameFilter, 
                      fileEndsWith=fileEndsWith,
                      excludeName=excludeName,
                      long=long, 
                      createdAfter=createdAfter)

def print_s3_files_flat(tree, nameFilter=None, fileEndsWith=None,
                        excludeName=None, long=False, createdAfter=None):
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
    long : bool
        Print all file info.
    createdAfter : datetime
        Only print files with creation data after.
    """
    # Print files at this level
    for f, info, date in zip(tree.get("files", []),
                             tree.get("info", []),
                             tree.get("date", [])):
        fname = f.split("/")[-1]
        if createdAfter is not None and date < createdAfter:
            continue
        if excludeName and excludeName in fname:
            continue
        if nameFilter and nameFilter not in fname:
            continue
        if fileEndsWith and not fname.endswith(fileEndsWith):
            continue
        if long:
            print(info, end=' ')
        print(f)

    # Recurse into subdirectories
    for dirname, subtree in tree.get("directories", {}).items():
        #if excludeName and excludeName in dirname:
        #   continue
        #if nameFilter and nameFilter not in dirname:
        #    continue
        print_s3_files_flat(subtree, nameFilter=nameFilter,
                            fileEndsWith=fileEndsWith,
                            excludeName=excludeName,
                            long=long,
                            createdAfter=createdAfter)

def run():
    
    initialLink, searchKeywords, printKeywords, flat = parseCommandLine()
    # Search for files
    s3files = list_s3_tree_recursive(initialLink, **searchKeywords)
    # Print the initial link, then print the rest
    print('\n', initialLink)
    if flat:
        print_s3_files_flat(s3files, **printKeywords)
    else:
        print_s3_tree(s3files, **printKeywords)

if __name__ == "__main__":
    run()