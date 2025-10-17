#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 14:35:25 2024

@author: ian
"""
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable


def autoScaleRange(myVar, percentile):
    '''
    If percentile less than 100, will select vmin as (100-percentile)
    and vmax as percentile of myVar, unless they fall out of the vmin and
    vmax bounds.
    Parameters
    ----------
    myVar : nparray
        Data being displayed.
    vmin : float
        Absolute minimum value.
    vmax : TYPE
        absolute maximum value.
    percentile : TYPE
        Clip data at (100-percentile) and percentile.

    Returns
    -------
    vmin, vmax - updated values based on percentiles.
    '''
    vmax = np.percentile(myVar[np.isfinite(myVar)], percentile)
    vmin = np.percentile(myVar[np.isfinite(myVar)], 100 - percentile)
    return vmin, vmax


def createDivider(ax, colorBarPosition='right', colorBarSize='5%',
                  colorBarPad=0.05):
    '''
    Create divider for color bar
    '''
    divider = make_axes_locatable(ax)  # Create space for colorbar
    return divider.append_axes(colorBarPosition, size=colorBarSize,
                               pad=colorBarPad)


def colorBar(pos, ax, colorBarLabel, colorBarPosition, colorBarSize,
             colorBarPad, labelFontSize, plotFontSize, extend='max',
             fontScale=1):
    '''
    Color bar for image
    '''
    # Create an divided axis for cb
    cbAx = createDivider(ax, colorBarPosition=colorBarPosition,
                         colorBarSize=colorBarSize,
                         colorBarPad=colorBarPad)
    # Select orientation
    orientation = {'right': 'vertical', 'left': 'vertical',
                   'top': 'horizontal',
                   'bottom': 'horizontal'}[colorBarPosition]

    cb = plt.colorbar(pos, cax=cbAx, orientation=orientation,
                      extend=extend)
    cb.set_label(colorBarLabel, size=labelFontSize * fontScale)
    cb.ax.tick_params(labelsize=plotFontSize * fontScale)
    if colorBarPosition in ['right', 'left']:
        cbAx.yaxis.set_ticks_position(colorBarPosition)
        cbAx.yaxis.set_label_position(colorBarPosition)
    elif colorBarPosition in ['top', 'tottom']:
        cbAx.xaxis.set_ticks_position(colorBarPosition)
        cbAx.xaxis.set_label_position(colorBarPosition)
