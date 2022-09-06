import numpy as np


from netCDF4 import Dataset

def dm_lwc(nw,lwc,rho):
    """
    returns dm[cm] as a function of nw[cm-4],
    lwc[g/m3] and rho_hydr[kg/m^3]
    """
    dm=(lwc*1e-3*4**4/(nw*np.pi*rho))**(0.25)
    return dm

from scipy.special import gamma as gam

def fmu(mu):
    return 6/4**4*(4+mu)**(mu+4)/gam(mu+4)
