import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import poisson, laplace, lognorm
from astropy.io import fits
from astropy.table import Table,join
import catalog_spect
import makemap
##colors
c6 = '#e41a1c'
c2 = '#377eb8'
c3 = '#4daf4a'
c4 = '#984ea3'
c5 = '#ff7f00'
c1 = '#66c2a5'
c7 = '#a65628'
c8 = '#f781bf'
c9 = '#999999'
colors = [c1,c2,c3,c4,c5,c6,c7,c8,c9]

#constants 
kpc = 3.08567758e21 #cm
Rsun = 8.3 #kpc
##Setting up some things for later
minflux,maxflux=-14,-6
bins = np.logspace(minflux,maxflux,25)
bins_mean = (bins[1:]+bins[:-1])/2.

#lb regions
lb = [(180,20),(180,5),(180,1),(20,20),(20,5),(20,1)]

##Set up bins for energy for fhl (10-1000GeV) and fgl (0.1-100GeV)
Efhl = np.logspace(1,3,100)#GeV
E3fgl_bins = np.logspace(-1,2,101)
E3fgl = (E3fgl_bins[1:]+E3fgl_bins[0:-1])/2.
E3fgl_dE = (E3fgl_bins[1:]-E3fgl_bins[0:-1])
mask_E3fgl = (E3fgl>=1)&(E3fgl<=100)
#bins in energy for our use (1-1000GeV)
Ebins = np.logspace(0,3,101)
E = (Ebins[1:]+Ebins[0:-1])/2.
dE = (Ebins[1:]-Ebins[0:-1])
mask_E3fhl = E>100
mask_E = (E>=1)&(E<=100)

###From tevcat 
Lmin_tev = 3.96e34
Lmax_tev = 1.79e38
#Function definitions
##Drawing random numbers
def rand_R(n,a,b,Rsun=Rsun):
    "Return random value of R, pdf is above, this is a gamma distribution"
    return np.random.gamma(a+1,1./b,n)*Rsun

def rand_R_norm(n,mu,sig):
    "return normal random"
    return np.random.normal(mu,sig,n)

def rand_R_SNI(n,a,Rsun=Rsun):
    return np.random.exponential(a,n)

def R_SNR_ferriere(R):
    sig1 = 4.8*np.exp(-(R-Rsun)/4.5)
    sig2 = 27*np.where(R < 3.7,3.55*np.exp(-((R-3.7)/2.1)**2),
                   np.exp(-(R**2-Rsun**2)/6.8**2))
    return sig1+sig2

###CDFs for R ferriere only
rspace = np.linspace(0,30,1000)
data_r = R_SNR_ferriere(rspace)
data_r = data_r/np.sum(data_r) #normalize
data_r_cdf = np.cumsum(data_r)#cdf

def ITS(N,x,cdf):
    """Function to do inverse transform sampling 
    
    Input:
    N: int, number of desired random samples
    x: array, possible returned random samples
    cdf: cdf function
    
    """
    smpl = np.zeros(N)
    u = np.random.rand(N)
    for i in range(N):
        #Inverse transform sampling, numerical
        ix= np.argmax(u[i]<=cdf)
        smpl[i]=x[ix]
    return smpl

def rand_z(n,H=0.33):
    """Return random z from laplace distribution=exp(-abs(z))
    n = number of randoms"""
    return np.random.laplace(0,H,n)

def z_SNR_ferriere(z):
    "z in pc"
    R1 = 7.3*np.exp(-np.abs(z)/325.)
    R2 = 50 * (0.79*np.exp(-(z/212.)**2)+0.21*np.exp(-(z/636.)**2))
    return R1+R2

###CDFs for z ferriere only
zspace = np.linspace(-3000,3000,1000)
data_z = z_SNR_ferriere(zspace)
data_z = data_z/np.sum(data_z) #normalize
data_z_cdf = np.cumsum(data_z)#cdf

def rand_z_SNR_ferriere(N):
    z_smpl = np.zeros(N)

    u = np.random.rand(N)
    for i in range(N):
        #Inverse transform sampling, numerical
        xx= np.argmax(u[i]<=data_z_cdf)
        z_smpl[i]=z[xx]
    return z_smpl/1000. #kpc

def rand_L(u,a,Lmin,Lmax):
    """Random L generator. Uses inverse transform sampling.
        u is a uniform random number between 0-1, a is the index of the
        power law distributiond desired.
        L = CDF^-1(u)
    """
    if a==1:
        print "a==1, error" #does not work for this..
        return 0
    return (u*(Lmax**(1-a)-Lmin**(1-a))+Lmin**(1-a))**(1/(1-a))

def norm_spec(spec,flux,dE,mask):
    """This function normalizes the spec(actually wrong_spec below) 
    to the integrated flux from 1-100GeV
    """
    wrong_spec = spec*flux
    norm = np.sum(wrong_spec[mask]*dE[mask])
    return wrong_spec*flux/norm

def read_catalog(fn,E,classes=['pwn','snr','spp']):
    """
    Reads the catalog, either 3FGL or 3FHL, and returns the average 
    spectrum of the desired sources.
    
    Returns flux,spec of desired source classes
    """
    if 'psch' in fn:
        fhl = True
        fgl = False
    else:
        fgl = True
        fhl = False
    data = fits.getdata(fn)
    tb = Table(data)
    if fgl:
        cls = np.asarray(tb['CLASS1'].tolist())
    elif fhl:
        cls = np.asarray(tb['CLASS'].tolist())
    cls = np.char.lower(np.char.rstrip(cls,' '))
    mask = np.zeros_like(cls,dtype=bool)
    lat = tb['GLAT']
    for m in classes:
        mask = (cls==m)|mask
        if 'lat' in m:
            mask = ((np.abs(lat)<5)&(cls==''))|mask
     
    if fgl:
        flux = tb[mask]['Flux1000']
    elif fhl:
        flux = tb[mask]['Flux']

    spec = catalog_spect.get_spec(fn,E) ###This function returns the spectra of the catalogue    
    spec = spec[mask]
    return flux,spec

def get_spec(Emin,Emax,index,flux,E):
    """
    For some flux between 1-100 GeV return power-law spectrum with index
    """

    if Emax < Emin:
        print "Emax < Emin"
        raise RuntimeError
    Efact = Emax**(1.-index)-Emin**(1-index)
    a = ((1-index)*(flux/Efact))
    return (a*(E[:,np.newaxis]**-index)).T

##This is the galaxy simulation
#Do the simulation of sources

def read_datfile(datFile):
    "Read dat file and return total gas spec (sum)"
    E_diff = datFile['emean']
    deltaE = datFile['delta_e']
    total_gas_spec = np.zeros_like(E_diff)
    total_gas_err = np.zeros_like(E_diff)
    mask_diff = E_diff > 1
    for d in datFile:
        if 'Gas' in d and not 'err' in d:
            total_gas_spec+=datFile[d]
            total_gas_err +=datFile[d+'_err']
    return total_gas_spec/deltaE,E_diff,total_gas_err/deltaE

pl_exp = lambda E, p1,p2,p3: p1*E**-p2*np.exp(p3*(1.-E))

def add_rand_spec(spec,flux,dE,mask):
    """
    spec is array with same length as flux and contains
    random spectra
    This functions normalizes the spectra to the correct value of flux
    """
    wrong_spec = flux[:,np.newaxis]*spec
    norm = np.sum(wrong_spec[:,mask]*dE[mask][np.newaxis,:],axis=1)
    good_spec = wrong_spec*(flux[:,np.newaxis]/norm[:,np.newaxis])
    return good_spec
                
def source_sim(N,a,Lmin,Lmax,specs,E,dE,thres=1e-9,spatialmodel='LorimerC',catalog='fgl'):
    #np.random.seed(0)
    """Sample from above distribution N times: r,z,L and calculate flux as measured at earth.
    New: add random spectrum to source 
    Return spec,flux,l,b
    """
#    np.random.seed(0)
    
    N = int(N) #no samples just to be sure 
    x0,y0,z0 = Rsun,0,0 #Location of us
    ##
    ##Drawing random numbers here
    u2 = np.random.rand(N)# For L
    phi_smpl = np.random.rand(N)*2*np.pi #uniform distribution for phi
    L_smpl = rand_L(u2,a,Lmin,Lmax)#random L
    
    if spatialmodel=='ferriere':
        R_smpl = ITS(N,rspace,data_r_cdf)
        z_smpl = ITS(N,zspace,data_z_cdf)/1000. #kpc
    elif spatialmodel=='LorimerC':
        alpha = 1.9
        beta = 5.0
        H = 0.18
        z_smpl = rand_z(N,H)
        R_smpl = rand_R(N,alpha,beta)
    elif spatialmodel=='LorimerS':
        alpha=0.2
        beta = 1.4
        H=0.33
        z_smpl = rand_z(N,H)
        R_smpl = rand_R(N,alpha,beta)
    elif spatialmodel== 'SNRGreen':
        alpha=1.09
        beta = 3.98
        H=0.083
        z_smpl = rand_z(N,H)
        R_smpl = rand_R(N,alpha,beta)
    else:
        print "spatialmodel not defined or wrong"
        raise RuntimeError
        return -1

    #distances from source to earth, calc flux
    x,y,z = R_smpl*np.cos(phi_smpl),R_smpl*np.sin(phi_smpl),z_smpl
    dist = np.sqrt(((x-x0))**2+((y-y0))**2+(z-z0)**2)
    distkpc = (dist*kpc)
    flux = L_smpl/(4*np.pi*distkpc**2)
    ##Calculate longitude and latitude
    l = np.arctan2((y0-y),(x0-x))
    b = np.arctan((z-z0)/np.sqrt((x-x0)**2+(y-y0)**2))
    l *= 180/np.pi
    b *= 180/np.pi
    
    #unresolved
    
#    #specs
    if catalog=='fgl':
        maskE = (E>=1)&(E<100)
    elif catalog=='fhl':
        maskE = (E>=10)
    
    randspec = specs[np.random.randint(0,specs.shape[0],N)]
    
    specs =  add_rand_spec(randspec,flux,dE,maskE)
    lb_arr = np.array((l,b))
    #specs =  get_spec(1,100,ind,flux,E)
    return specs,flux,lb_arr
