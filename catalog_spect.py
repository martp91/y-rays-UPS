import numpy as np
import pylab as plt
import optparse
from astropy.io import fits
from astropy.table import Table

def get_spec(filename):
    if 'psch' in filename:
        fhl = True
        fgl = False
    else:
        fgl = True
        fhl = False
    try:
        data = fits.getdata(filename)
    except:
        print "does file exist?"
    tb = Table(data)
    if fgl:
        E = np.logspace(-1,2,100)#GeV
    if fhl:
        E = np.logspace(1,3,100)#GeV
    ST = np.char.rstrip(tb['SpectrumType'],' ')
    spect = np.zeros((len(ST),len(E)))
    #mask_pwl = (ST=='PowerLaw')
    #mask_log = (ST=='LogParabola')
    #mask_Exp = (ST=='PLExpCutoff')
    #mask_SuperExp = (ST=='PLSuperExpCutoff')
    #list_mask = [mask_pwl,mask_log,mask_Exp]
    #list_p = ['Flux_Density','Pivot_Energy','PowerLaw_Index','Spectral_Index','beta','Cutoff','Exp_Index']
    func_pwl = lambda E,p: p[0]*(E/p[1])**-p[2]
    func_log = lambda E,p: p[0]*(E/p[1])**(-p[3]-p[4]*np.log(E/p[1]))
    func_Exp = lambda E,p: p[0]*(E/p[1])**(-p[3])*np.exp((p[1]/p[5])**p[6]-(E/p[5])**p[6])
    list_func=[func_pwl,func_log,func_Exp]
    if fhl:
        fd = tb['Flux_Density']
        pe = tb['Pivot_Energy']
        pi = tb['PowerLaw_Index']
        si = tb['Spectral_Index']
        beta = tb['beta']
    if fgl:
        fd = tb['Flux_Density']*1000. #GeV
        pe = tb['Pivot_Energy']/1000. #GeV
        pi = tb['PowerLaw_Index']
        si = tb['Spectral_Index']
        beta = tb['beta']
        co = tb['Cutoff']/1000.
        ei = tb['Exp_Index']

    #for i,m in enumerate(list_mask):
    #    spect[m,:]= list_func[i](E,[fd[m],pe[m],pi[m],si[m],beta[m],co[m],ei[m]])
    #[fd[i],pe[i],pi[i],si[i],beta[i],co[i],ei[i]

    for i,entry in enumerate(tb):
        if ST[i]=='PowerLaw':
            spect[i,:]=func_pwl(E,[fd[i],pe[i],pi[i]])
        elif ST[i]=='LogParabola':
            spect[i,:]=func_log(E,[fd[i],pe[i],pi[i],si[i],beta[i]])
        elif ST[i]=='PLExpCutoff':
            spect[i,:]=func_Exp(E,[fd[i],pe[i],pi[i],si[i],beta[i],co[i],ei[i]])
        elif ST[i]=='PLSuperExpCutoff':
            spect[i,:]=func_Exp(E,[fd[i],pe[i],pi[i],si[i],beta[i],co[i],ei[i]])
        else:
            print "missed:", ST[i] 
        #plt.loglog(E,spect[i,:])
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



