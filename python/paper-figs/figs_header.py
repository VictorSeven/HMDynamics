#!/usr/bin/python

from matplotlib import rcParams

TICKFONTSIZE = 10
LABELFONTSIZE = 12

def mm_to_inch(mm):
    """
    mm: Lenght in millimeters to be converted to inches
    """
    return mm / 25.4

def figsize(width_mm = 190, ratio = 1.618):
    """
    width_mm: figure width in millimeters, default to 190mm (Elsevier single column)
    ratio: width/height relationship  
    """
    w = mm_to_inch(width_mm)
    h = w / ratio
    return (w,h)

def set_up_figure():
    """
    Set up font family and size for figures
    """

    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = "Helvetica"

    rcParams["axes.labelsize"] = LABELFONTSIZE
    rcParams["legend.fontsize"] = LABELFONTSIZE
    rcParams["xtick.labelsize"] = TICKFONTSIZE
    rcParams["ytick.labelsize"] = TICKFONTSIZE
