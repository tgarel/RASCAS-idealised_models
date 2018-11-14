def write_fits_image(image,pixstep,fitsFileName):
    import astropy.io.fits as fits
    import astropy.wcs as pywcs
    from datetime import datetime

    primary_header=fits.Header()
    data_header=fits.Header()
    wcs = pywcs.WCS(naxis=2)
    wcs.wcs.crpix = np.array([(image.shape[0]+1)/2.,(image.shape[1]+1)/2.])
    wcs.wcs.crval = np.array([0.,0.])
    wcs.wcs.ctype = ['x', 'y']
    wcs.wcs.cunit = ['arcsec', 'arcsec']
    wcs.wcs.crota=[0.0,0.0]
    wcs.wcs.cdelt=[pixstep,pixstep]
    hdrwcs = wcs.to_header()
    primary_header['date'] = (str(datetime.now()), 'creation date')
    primary_header['author'] = ('J.Blaizot, L. Michel-Dansac', 'SPHINX simulation')
    hdulist = fits.HDUList([fits.PrimaryHDU(header=primary_header)])
    keys = set(data_header.keys()) - set(hdrwcs.keys())
    for card in data_header.cards:
        if card.keyword not in keys:
            continue
        hdrwcs[card.keyword] = (card.value, card.comment)
    hdrwcs['FUNITS']='erg/s/cm2/A/arcsec2'
    hdulist.append(fits.ImageHDU(name='DATA',data=image,header=hdrwcs))
    hdulist.writeto(fitsFileName, clobber=True,output_verify='silentfix')

    
###########33

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
plt.ion()
import rascasRun as rr
import lya_utils as lya

surveyName = ['F115W', 'F150W', 'F200W', 'F277W', 'F356W', 'F444W']

rascasDir      ='/cral2/blaizot/RASCAS_test/halo41361/'
DomDumpDir     = 'CDD_HI_dust'
ramsesDir      =''
ramsesTimestep =1



# check spectra 
if False:
    for f in surveyName:
        surv = rr.RascasSurvey(f,rascasDir,DomDumpDir,ramsesDir,ramsesTimestep)
        surv.load(savePhots=True)
        x1,s1 = surv.pesc.spectrum(frame='obs')
        x2,s2 = surv.p.spectrum(frame='ic')
        plt.plot(x1,s1)
        plt.plot(x2,s2)
    
# compute images ...
from minirats.HaloFinder.py import haloCatalog as hC

ramsesDir      = '/cral2/sphinx/05_F1000/02_IC20_BP/'
ramsesTimestep = 183 
hcat = hC.haloCatalog(ramsesDir,ramsesTimestep,zoom=False,HaloDir="%s/Halos/"%(ramsesDir))
redshift = hcat.info['redshift']
from astropy.cosmology import FlatLambdaCDM, z_at_value
import astropy.units as u
cosmo = FlatLambdaCDM(H0=hcat.info['H0'],Om0=hcat.info['omega_m'])
lumDist = cosmo.luminosity_distance(hcat.info['redshift'])
lumDist_cm = lumDist.cgs.value
a =  cosmo.arcsec_per_kpc_proper(hcat.info['redshift'])
pc2cm = 3.086e18
cm2pc = 1./pc2cm
unit_l_pkpc = hcat.info['unit_l'] * cm2pc / 1000. # converts comoving code units to phys. kpc
unit_l_arcsec = unit_l_pkpc * a.value


class snap_class(object):
    def __init__(self,redshift,unit_l_arcsec,lumdist_cm):
        self.redshift      = redshift
        self.unit_l_arcsec = unit_l_arcsec
        self.lumDist_cm    = lumdist_cm

snap = snap_class(redshift,unit_l_arcsec,lumDist_cm)


from nircam import nircamFilter
import numpy as np
for f in surveyName:
    surv = rr.RascasSurvey(f,rascasDir,DomDumpDir,ramsesDir,ramsesTimestep)
    surv.load()
    lbdapix = 20. #[A]
    filterCurve = nircamFilter(f)

    lbdamin_restframe=filterCurve.lambda_min/(1+redshift)
    lbdamax_restframe=filterCurve.lambda_max/(1+redshift)

    # NIRCAM parameters :
    pixsize = 0.01 # [arcsec]
    nxybins = 500

    for k in [[1.,0,0],[0,1,0],[0,0,1]]:

    
        cube,lmin,nlbins,xmin,ymin,nxybins = surv.make_mock_cube(k,snap,thetamax=20.,nxybins=nxybins,xc=None,yc=None,zc=None,\
                                                                    pixsize=pixsize,lbdapix=lbdapix,\
                                                                    lbdamin_restframe=lbdamin_restframe,lbdamax_restframe=lbdamax_restframe)
        lbda = np.linspace(lmin,lmin+nlbins*lbdapix,nlbins)
    
        throughput  = np.interp(lbda,filterCurve.wl,filterCurve.throughput)  #Interpolate to common wavelength axis
        throughput = throughput * lbda
        
        im = np.average(cube,axis=0,weights=throughput,returned=False)  # erg / s / cm2 / arcsec2 / A
        # x lambda^2/c pour avoir Fnu -> -2.5log10(Fnu)-48.6    
        
        write_fits_image(im,pixsize,'%s/%s/%.2f_%.2f_%.2f.fits'%(rascasDir,f,k[0],k[1],k[2]))

        if False:
            plt.figure()
            plt.imshow(im.T,interpolation='nearest',norm=LogNorm(),vmin=1e-23,vmax=1e-18,
                        origin='lower',cmap='gist_gray_r')
            plt.xlabel('arcsec')
            plt.colorbar()







    
    
    


    
