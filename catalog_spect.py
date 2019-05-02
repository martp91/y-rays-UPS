import numpy as np
import pylab as plt
import optparse
from astropy.io import fits
from astropy.table import Table

def get_spec(filename,E):
    if 'psch' in filename:
        fhl = True
        fgl = False
    else:
        fgl = True
        fhl = False
    fgl4 = False
    if '8year' in filename or 'v18' in filename:
        fgl4 = True
    try:
        data = fits.getdata(filename)
    except:
        print "does file exist?"
    tb = Table(data)
    ST = np.char.rstrip(tb['SpectrumType'],' ')
    spect = np.zeros((len(ST),len(E)))
    
    func_pwl = lambda E,p: p[0]*(E/p[1])**-p[2]
    func_log = lambda E,p: p[0]*(E/p[1])**(-p[2]-p[3]*np.log(E/p[1]))
    func_Exp = lambda E,p: p[0]*(E/p[1])**(-p[3])*np.exp((p[1]/p[5])**p[6]-(E/p[5])**p[6])

    if fgl4:
        func_subExp = lambda E,p: p[0]*(E/p[1])**-p[2]*np.exp(p[3]*(p[1]**p[4]-E**p[4]))

    list_func=[func_pwl,func_log,func_Exp]
    if fhl:
        fd = tb['Flux_Density']
        pe = tb['Pivot_Energy']
        pi = tb['PowerLaw_Index']
        si = tb['Spectral_Index']
        beta = tb['beta']
    if fgl and not fgl4:
        fd = tb['Flux_Density']*1000 #GeV
        pe = tb['Pivot_Energy']/1000. 
        pi = tb['PowerLaw_Index']
        si = tb['Spectral_Index']
        beta = tb['beta']
        co = tb['Cutoff']/1000.
        ei = tb['Exp_Index']
    if fgl4:
        fd = tb['PL_Flux_Density']*1000.
        pe = tb['Pivot_Energy']/1000.
        pi = tb['PL_Index']
        li = tb['LP_Index'] #logparabola index
        lpbeta = tb['LP_beta']
        plec_i = tb['PLEC_Index'] #low E expcutoff
        plec_ef = tb['PLEC_Expfactor'] #expfactor a
        plec_ei = tb['PLEC_Exp_Index'] #expindex b


    for i,entry in enumerate(tb):
        if ST[i]=='PowerLaw':
            spect[i,:]=func_pwl(E,[fd[i],pe[i],pi[i]])
        elif ST[i]=='LogParabola' and not fgl4:
            spect[i,:]=func_log(E,[fd[i],pe[i],si[i],beta[i]])
        elif ST[i] =='LogParabola' and fgl4:
            spect[i,:] = func_log(E,[fd[i],pe[i],li[i],lpbeta[i]])
        elif ST[i]=='PLExpCutoff':
            spect[i,:]=func_Exp(E,[fd[i],pe[i],pi[i],si[i],beta[i],co[i],ei[i]])
        elif ST[i]=='PLSuperExpCutoff' and not fgl4:
            spect[i,:]=func_Exp(E,[fd[i],pe[i],pi[i],si[i],beta[i],co[i],ei[i]])
        elif ST[i]=='PLSuperExpCutoff2' and fgl4:
            spect[i,:] = func_subExp(E,[fd[i],pe[i],plec_i[i],plec_ef[i],plec_ei[i]])
        else:
            print "missed:", ST[i] 
        #plt.loglog(E,spect[i,:])
    #plt.show()
    return spect

#TODO make E variable???    
def average_spec(spec):
    from scipy.optimize import curve_fit
    "Return ave spec and best fit pwl index"
    _E = np.logspace(2,5,100)#MeV

    fit_func = lambda E,a,b: a-b*np.log10(E/1000.)
    ave_spec = np.average(spec,axis=0)
    fit,cov = curve_fit(fit_func,_E,np.log10(ave_spec),p0=[-12,2.2])

    return ave_spec,fit[1] 



