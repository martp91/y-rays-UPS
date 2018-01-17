import ast
import numpy as np

def open_dat(fname):
    """open dat file return dictionary with spectra"""

    open_file = open(fname, 'r')
    comps = ['ICS 0-3kpc', 'ICS 3-8.3kpc','ICS 8.3-50kpc','ICSinner','ICSouter','emean','delta_e','Gas ring I','Gas ring II','Gas ring III', 'ICS', 'Bubbles', 'data', 'IGRB','PS','511 keV temp','Gas ring 1', 'Gas ring 2','Gas ring 3', 'Gas ring 4', 'Gas ring 5', 'Gas ring 6', 'Gas ring 7', 'Gas ring 8', 'Gas ring 9', 'Gas ring 25'] #Change according to components
    err = [c+'_err' for c in comps]
    comps = comps+err
    d = {}
    for line in open_file:
        entry = line.split(':')
        if entry[0] in comps:
            d[entry[0]]= np.array(ast.literal_eval(entry[1].strip()))
    return d
