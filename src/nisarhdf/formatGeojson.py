#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 08:46:28 2024

@author: ian
"""


def myIndent(n, x, formattedString):
    for space in [' '] * n * 2:
        formattedString += space
    return formattedString


def formatGeojson(myString):
    '''
    Format a geojson string for readability.

    Parameters
    ----------
    myString : str
        GeoJson unformatted string

    Returns
    -------
    formatted geoJson string

    '''
    # indent function
    #
    # Start bracket counds at 0
    braceCount, bracketCount = 0, 0
    # Remote space
    myString = myString.replace(' \"', '\"')
    #
    # Start with empty string
    formattedString = ''
    # And add characters back with additional formatting as neeed.
    for x in myString:
        # Add return and indent if not in list (leftBracket==0)
        if x in ',' and bracketCount == 0:
            formattedString += x
            formattedString += '\n'
            formattedString = myIndent(braceCount, x, formattedString)
        # Update bracket count if left bracket
        elif x == '[':
            bracketCount += 1
            formattedString += x
        # decrement left bracket count
        elif x == ']':
            bracketCount -= 1
            formattedString += x
        # Update brace count if left brace and add return and indent
        elif x == '{':
            braceCount += 1
            formattedString += x
            formattedString += '\n'
            formattedString = myIndent(braceCount, x, formattedString)
        # Decrement braceCount and add new line and indent
        elif x == '}':
            formattedString += '\n'
            braceCount -= 1
            formattedString = myIndent(braceCount, x, formattedString)
            formattedString += x
        # Regulary character so append
        else:
            formattedString += x
    return formattedString
