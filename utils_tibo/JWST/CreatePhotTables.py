import argparse
import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from scipy.io import FortranFile
import numpy as np
import os
import sys

# RASCAS 
import lya_utils as lya


def main(args):

    # Read SSP spectra from sedDir 
    ssp = readRamsesSEDs(args.sedDir)
    # generate tables
    if args.specType == 'Table':
        gen_tabulated_spectra(args.llow,args.lup,ssp,args.outputfile,dlambda=args.dlambda,comment=args.comment)
    elif args.specType == 'Mono':
        gen_monochromatic_tables(args.lambda0,ssp,args.outputfile,comment=args.comment)
    elif args.specType == 'Lya':
        gen_LyaLum_tables(ssp,args.outputfile,comment=args.comment)
    elif args.specType == 'PowLaw':
        gen_PowLaw_tables(ssp,args.llow,args.lup,args.lambda0,args.outputfile,comment=args.comment,plotFile=args.plotFile,removeAbsLines=args.removeAbsLines)
        

def gen_tabulated_spectra(llow,lup,ssp,outputfile,dlambda=0.01,comment=''):
    # inputs:
    # - llow : min wavelength [A]
    # - lup  : max wavelength [A]
    # - ssp  : SSP spectra (as read by readRamsesSEDs)
    # - outputfile = name of table file
    # - dlambda : binning of tables [A] (0.01 is a good choice)
    # - comment : an optional sentence to add in first line of output file (e.g. giving library source).

    # look up for llow and lup in ssp's lambdas 
    ilow = np.where(abs(ssp['lambdaBins']-llow) == min(abs(ssp['lambdaBins']-llow)))[0][0]-1
    llow = ssp['lambdaBins'][ilow]
    iup  = np.where(abs(ssp['lambdaBins']-lup) == min(abs(ssp['lambdaBins']-lup)))[0][0]+1
    lup = ssp['lambdaBins'][iup]
    flux = ssp['spectra'][ilow:iup,:,:] * ssp['Lsun_cgs'] # [erg/s/Msun/A]
    lbda = ssp['lambdaBins'][ilow:iup]  # [A]

    # binning 
    ngood = int((lup - llow)/dlambda)
    nage  = len(ssp['ageBins'])
    nmet  = len(ssp['ZBins'])

    f = open(outputfile,'w')
    f.write("# Tabulated SEDs, from %f to %f A. %s \n"%(llow,lup,comment))
    f.write("# First line is nages and nmets. Then come ages [Myr], then metallicities [absolute] (lines 2 and 3). \n")
    f.write("# Line 4 is nbins. Following lines give nb of photons [#/s], lambda(x) with x regularly binned (nbins) in [0,1] (one line per age and metallicity). \n")
    f.write("%i %i \n"%(nage,nmet))
    for iage in range(nage): f.write("%.8e "%(ssp['ageBins'][iage]*1e-6))
    f.write("\n")
    for imet in range(nmet): f.write("%.8e "%(ssp['ZBins'][imet]))
    f.write("\n")
    f.write("%i \n"%(ngood))

    # compute normalisation of each SED (nb of photons emitted in wavelength range)
    for imet in range(nmet):
        for iage in range(nage):
            # get N(lambda)=F_lambda / (h nu)
            fl    = np.ravel(flux[:,iage,imet])  # [erg/s/A/Msun]
            fl    = fl * lbda * 1.0e-8 / lya.h_cgs / (lya.clight) # [#/s/A/Msun]
            
            # rebin specs before integration
            lll = np.linspace(lbda[0],lbda[-1],ngood)
            fff = np.interp(lll,lbda,fl)
            
            # compute the cumulative proba distribution of photons... (P(<lambda))
            p = [0]
            l = [lbda[0]]
            pp = 0.
            for i in np.arange(1,len(lll)):
                pp = pp + 0.5*(fff[i-1]+fff[i])*(lll[i]-lll[i-1])
                p.append(pp)
                l.append(lll[i])
    
            f.write("%.8e "%(pp)) # total nb of photons [#/s/Msun] 
        
            x = np.ravel(l)
            y = np.ravel(p/max(p))   # P(<lambda)
            # rebin to get the reciprocal in regular bins
            ygood = np.linspace(0,1,ngood)
            xgood = np.interp(ygood,y,x)
            for xg in xgood: f.write("%.8e "%(xg))
            f.write("\n")
        
    f.close()


def gen_monochromatic_tables(lambda0,ssp,outputfile,comment=''):
    # inputs :
    # - lambda0 : wavelength at which we want the flux [A]
    # - ssp  : SSP spectra (as read by readRamsesSEDs)
    # - outputfile = name of table file
    # - comment : an optional sentence to add in first line of output file (e.g. giving library source).

    i0 = np.where(abs(ssp['lambdaBins']-lambda0) == min(abs(ssp['lambdaBins']-lambda0)))[0][0]
    lambda0 = ssp['lambdaBins'][i0]
    f0  = ssp['spectra'][i0,:,:] * ssp['Lsun_cgs']  # erg/s/A/Msun
    nu0 = lya.clight / (lambda0*1e-8) # [Hz]
    n0  = f0 / lya.h_cgs / nu0  # nb of photons / s / A / Msun 
    
    # save these to some ascii file 
    f = open(outputfile,'w')
    f.write("# Number of photons emitted /s/Msun/A, at %f A. %s \n"%(lambda0,comment))
    f.write("# First line is nages and nmets. Then come ages [Myr] (line 2), metallicities [absolute] (line 3), and wavelength [A] (line 4) \n")
    f.write("# Then each line has N_%i [#/s/Msun/A] for all ages and one metallicity \n"%(lambda0))
    nage = len(ssp['ageBins'])
    nmet = len(ssp['ZBins'])
    f.write("%i %i \n"%(nage,nmet))
    for iage in range(nage): f.write("%.8e "%(ssp['ageBins'][iage]*1e-6))
    f.write("\n")
    for imet in range(nmet): f.write("%.8e "%(ssp['ZBins'][imet]))
    f.write("\n")
    f.write("%.8e \n"%(lambda0))
    for imet in range(nmet):
        for iage in range(nage):f.write("%.8e "%(n0[iage,imet]))
        f.write("\n")
    f.close()

    
def gen_LyaLum_tables(ssp,outputfile,comment=''):
    # inputs :
    # - ssp  : SSP spectra (as read by readRamsesSEDs)
    # - outputfile = name of table file
    # - comment : an optional sentence to add in first line of output file (e.g. giving library source).

    f = open(outputfile,'w')
    f.write("# Number of Lya photons emitted /s/Msun. %s \n"%(comment))
    f.write("# First line is nages and nmets. Then come ages [Myr] (line 2), metallicities [absolute] (line 3), and central wavelength [A] (line 4) \n")  
    f.write("# Then each line has N_Lya [#/s/Msun] for all ages and one metallicity \n")
    nage = len(ssp['ageBins'])
    nmet = len(ssp['ZBins'])
    f.write("%i %i \n"%(nage,nmet))
    for iage in range(nage): f.write("%.8e "%(ssp['ageBins'][iage]*1e-6))
    f.write("\n")
    for imet in range(nmet): f.write("%.8e "%(ssp['ZBins'][imet]))
    f.write("\n")
    f.write("%.8e \n"%(lya.lambda0*1e8))  # line-center wavelength [A]
    # comute number of Lyman continuum photons
    lylim  = 912. # [Angstrom]
    ilylim = np.where(abs(ssp['lambdaBins']-lylim) == min(abs(ssp['lambdaBins']-lylim)))[0][0]
    lbds   = ssp['lambdaBins'][:ilylim+1]
    for imet in range(nmet):
        x,y = [],[]
        for iage in range(nage):
            fly  = np.ravel(ssp['spectra'][0:ilylim+1,iage,imet]) * ssp['Lsun_cgs'] # [erg/s/A/Msun]
            fly  = fly * lbds
            nlyc = np.trapz(fly,lbds)
            nlyc = nlyc * 1.0e-8 / lya.h_cgs / (lya.clight) # [#/s/Msun]
            nlya = nlyc * 0.68  # nb of Lya photons
            x.append(ssp['ageBins'][iage])
            y.append(np.log10(nlya))
            f.write("%.8e "%(nlya))
        f.write("\n")
    f.close()


def gen_PowLaw_tables(ssp,llow,lup,lambda0,outputfile,comment='',plotFile=None,removeAbsLines=False):
    # inputs :
    # - ssp  : SSP spectra (as read by readRamsesSEDs)
    # - llow : min wavelength [A]
    # - lup  : max wavelength [A]
    # - lambda0 : pivot wavelength [A]
    # - outputfile = name of table file
    # - comment : an optional sentence to add in first line of output file (e.g. giving library source).
    # - plotFile : optional filename for a plot
    # removeAbsLines : use some iterative scheme to filter out strong abs. lines. 
    
    # look up for llow and lup in ssp's lambdas 
    ilow = np.where(abs(ssp['lambdaBins']-llow) == min(abs(ssp['lambdaBins']-llow)))[0][0]-1
    llow = ssp['lambdaBins'][ilow]
    iup  = np.where(abs(ssp['lambdaBins']-lup) == min(abs(ssp['lambdaBins']-lup)))[0][0]+1
    lup = ssp['lambdaBins'][iup]
    flux = ssp['spectra'][ilow:iup,:,:] * ssp['Lsun_cgs'] # [erg/s/Msun/A]
    loglambda = np.log10((ssp['lambdaBins'])[ilow:iup]) # [A]
    logf      = np.log10(flux)

    
    # save these to some ascii file 
    f = open(outputfile,'w')
    f.write("# Continuum power-law fits, from %f to %f A, in the form F_l = F_0 * (l/l0)^beta, with l0 = %f. %s\n"%(llow,lup,lambda0,comment))
    f.write("# First line is nages and nmets. Then come ages [Myr], then metallicities [absolute] (lines 2 and 3), and l0 [A] (line 4) \n")
    f.write("# Then consecutive lines (corresponding to a single metallicity) have F_0 [erg/s/A/Msun] for all ages, beta for all ages. \n")
    nage = len(ssp['ageBins'])
    nmet = len(ssp['ZBins'])
    f.write("%i %i \n"%(nage,nmet))
    for iage in range(nage): f.write("%.8e "%(ssp['ageBins'][iage]*1e-6))
    f.write("\n")
    for imet in range(nmet): f.write("%.8e "%(ssp['ZBins'][imet]))
    f.write("\n")
    f.write("%.8e \n"%(lambda0))

    if plotFile is not None:
        # prepare plot 
        plt.figure(figsize=(12,10))
        import matplotlib.cm as cm
        colors = cm.rainbow(np.linspace(0, 1, 221))

    for imet in range(nmet):
        x,y = [],[]
        F_0s = ""
        betas = ""
        for iage in range(nage):
            p = np.polyfit(loglambda,logf[:,iage,imet],deg=1)
            fit = 10.**( p[1] + p[0]*loglambda) 
            if removeAbsLines:
                # remove the strong Lya absorption in the SEDs.
                ii = np.where(flux[:,iage,imet] > 0.8*fit)[0]
                p = np.polyfit(loglambda[ii],logf[ii,iage,imet],deg=1)
                fit = 10.**( p[1] + p[0]*loglambda) 
                # again, with stronger constraint
                ii = np.where(flux[:,iage,imet] > 0.95*fit)[0]
                p = np.polyfit(loglambda[ii],logf[ii,iage,imet],deg=1)
                fit = 10.**( p[1] + p[0]*loglambda) 
            # extract fit parameters
            beta   = p[0]
            logF_0 = p[1]+beta*np.log10(lambda0)
            F_0    = 10.**logF_0  # in erg/s/A/Msun
            F_0s = "%s %.8e"%(F_0s,F_0)
            betas = "%s %.8e"%(betas,beta)

            if plotFile is not None:
                if imet == 1 and iage%10==0:
                    plt.plot(ssp['lambdaBins'],ssp['spectra'][:,iage,imet]*ssp['Lsun_cgs'],color=colors[iage])
                    fit = F_0 * (10.**loglambda/lambda0)**beta
                    plt.plot(10.**loglambda,fit,color=colors[iage],linewidth=3,
                        alpha=0.6,label="%.1e Myr"%(ssp['ageBins'][iage]*1e-6))
            
        f.write("%s \n"%(F_0s))
        f.write("%s \n"%(betas))
        
    f.close()
    if plotFile is not None:
        plt.xlim(llow,lup)
        plt.ylim(4e26,4e33)
        plt.xlabel(r'$\lambda \ [\AA]$',fontsize=15)
        plt.ylabel(r'$F_{\lambda} \ [erg / s / \AA / M_\odot]$',fontsize=15)
        plt.yscale('log')
        plt.axvline(1216,linestyle='--',alpha=0.3,color='red',linewidth=4,label=r'$\lambda_{Ly\alpha}$')
        plt.legend(loc='upper right')
        plt.savefig(plotFile)
        plt.close()
        
    
# Function taken fron ramses utils, by Trebitsch & Rosdahl.
def readRamsesSEDs(sedDir):
    """Read SED in ramses format and return
    Parameters:
    ----------------------------------------------------------------------
    sedDir: Directory containing the SED tables
    """
    # Read metallicity bins
    ZFile = open(sedDir+'/metallicity_bins.dat', 'r')
    nZ = eval(ZFile.readline())
    ZBins = []
    for Z in range(0,nZ): ZBins.append(eval(ZFile.readline()))
    ZFile.close()

    # Read age bins
    ageFile = open(sedDir+'/age_bins.dat', 'r')
    nAge = eval(ageFile.readline())
    ageBins = []
    for age in range(0,nAge): ageBins.append(eval(ageFile.readline()))
    ageFile.close()

    # Read wavelength bins and spectra
    # spectra are in [Lsun / A / Msun]
    # with Lsun = 3.826e33 erg/s and Msun = 2e33g
    sedFile = FortranFile(sedDir+'/all_seds.dat','r')
    nLambda = sedFile.read_ints()[0]
    lambdaBins = sedFile.read_reals()
    spectra = np.empty([nLambda,nAge,nZ])
    for iZ in range(0,nZ):
        for iAge in range(0,nAge):
            spectrum = sedFile.read_reals()
            spectra[:,iAge,iZ] = spectrum  

    Lsun_cgs = 3.826e33
    Msun_cgs = 2e33
    return {'ZBins':ZBins, 'ageBins':ageBins, 'lambdaBins':lambdaBins,'spectra':spectra,'Lsun_cgs':Lsun_cgs,'Msun_cgs':Msun_cgs}



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-sed","--sedDir", help="path to SSP spectra ",default="./")
    parser.add_argument("-specType","--specType",help="can be Mono, PowLaw, Lya, Table (default)",default="Table")
    parser.add_argument("-llow","--llow",help="min wavelength [A] (for PowLaw or Table specType)",default=1100)
    parser.add_argument("-lup","--lup",help="max wavelength [A] (for PowLaw or Table specType)",default=1300)
    parser.add_argument("-outputfile","--outputfile",help="output file name",default="table.txt")
    parser.add_argument("-comment","--comment",help="comment to append to first line of output file",default="")
    parser.add_argument("-dlamnda","--dlambda",help="bin size [A] for SSP spectra tables (default 0.01)",default=0.01)
    parser.add_argument("-lambda0","--lambda0",help="wavelength [A] for specType Mono or PowLaw",default=1200)
    parser.add_argument("-plotFile","--plotFile",help="optional filename for plot with specType == PowLaw",default=None)
    parser.add_argument("-removeAbsLines","--removeAbsLines",help="filter out absorption lines if specType == PowLaw",action="store_true")
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    main(args)
    
