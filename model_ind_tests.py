def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn
import sys as sys

from gapp import gp, dgp, covariance
import pickle
import numpy as np
from numpy import array,concatenate,loadtxt,savetxt,zeros
import matplotlib.pyplot as plt
from matplotlib import rc

# ===================================================================================================================
# Code for performing model-independent tests using Euclid and SKA simulated H(z) measurements from radial BAO mode
# Written by Carlos Bengaly around 2019-2020, adapted in July 2021 
# Contact: carlosap87@gmail.com or carlosbengaly@on.br; github/astrobengaly
# ===================================================================================================================

# ================================================================================================================================================================= 
# Receiving the number of data points and name of the survey containing Hz measurements as inputs
# data-sets currently supported: 
# - real Hz measurements from galaxy age cosmic chronometers (CC) and radial BAO mode
# - Euclid and SKA simulations for radial BAO mode
# For real data: nz1 refers to the number of CC data points, nz2 to the BAO ones. Put 0 in nz2 in case you just want the former 
# For simulations: nz1 refers to the number of Euclid or SKA band 1 data points, nz2 to the SKA band 2. Put 0 in nz2 if you do not want to include band 2 simulations
# survey refers to 'euclid', 'ska', or 'cc'
# =================================================================================================================================================================

nz1=int(sys.argv[1])
nz2=int(sys.argv[2])
survey=sys.argv[3]

if __name__=="__main__":

    ##looping through nsim - not necessary for the time being
    #for n in range(nsim):
    
        # ====== Loading data
        filename = survey+'_'+str(nz1)+'+'+str(nz2)+'pts'
        (Z,Hz,Sigma,hid) = loadtxt('hz_'+filename+'.dat',unpack='True')

        # ======= Setting the redshift range where the reconstructions will be built
        zmin = 0.
        zmax = np.max(Z)

        #======== Gaussian process: reconstructing H(z) from loaded data, starting from zmin = 0 all the way to zmax using 100 bins
        g = gp.GaussianProcess(Z,Hz,Sigma,covfunction=covariance.SquaredExponential,cXstar=(zmin,zmax,500))
        dg = dgp.DGaussianProcess(Z,Hz,Sigma,covfunction=covariance.SquaredExponential,cXstar=(zmin,zmax,500))
        (rec,theta) = g.gp(thetatrain='False')
        (drec,dtheta) = dg.dgp(thetatrain='False')
        (d2rec,theta) = dg.d2gp(thetatrain='False')
                        
        #======== calculating covariances between hz, h'z and h''z at points Zstar.
        fcov = dg.f_covariances(fclist=[0,1])
        #savetext = ('cov_'+filename+'.dat', fcov)
        #f = open('cov_'+filename+'.dmp','wb')
        #pickle.dump(fcov,f)
        #f.close()
        
        #========= associating the reconstructed quantities with named variables
        H0 = 67.36
        n_start = 0
        zrec = rec[n_start:,0]
        hzrec = rec[n_start:,1]
        ezrec = rec[n_start:,1]/H0
        sighzrec = rec[n_start:,2]      
        sigezrec = rec[n_start:,2]/H0     
        dhzrec = drec[n_start:,1]
        dezrec = drec[n_start:,1]/H0
        sigdhzrec = drec[n_start:,2]
        sigdezrec = drec[n_start:,2]/H0
        sighdhzrec = fcov[n_start:,:,]
                
        #======== calculating the deceleration parameter qz and its uncertainty
        qzrec = ( (1.+zrec)*dhzrec/hzrec )-1.
        sigqzrec = np.sqrt( (sighzrec/hzrec)**2. + (sigdhzrec/dhzrec)**2. - 2.*sighdhzrec[:,1,0]/(hzrec*dhzrec) )*(1.+qzrec)
        
        # ========== omz and Lmz reconstructed from hz: flat LCDM null test where om = wm and Lm = 0 if the underlying Cosmology is flat LCDM
        om = (ezrec**2.-1.)/((1.+zrec)**3.-1.)
        sigom = (2.*ezrec*sigezrec)/((1.+zrec)**3.-1.) 
        Lm = 3.*((1.+zrec)**2.)*(1.-ezrec**2.) + 2.*zrec*(3 + 3.*zrec + zrec**2.)*ezrec*dezrec
        sigLm = np.sqrt( ( 3.*((1.+zrec)**2.)*(2.*ezrec) + 2.*zrec*(3 + 3.*zrec + zrec**2.)*dezrec )**2.*sigezrec**2. + ( (2.*zrec*(3 + 3.*zrec + zrec**2.)*ezrec )**2.*sigdezrec**2. ) )
        
        #====== printing the H0 and q0 values
        print 'z=', zrec[0], ' H0=', hzrec[0], ' sigH0=',  sighzrec[0], ' sigH0/H0 (%)=', (sighzrec[0]/hzrec[0])*100.
        print 'z=', zrec[0], ' q0=', qzrec[0], ' sigq0=',  sigqzrec[0], ' sigq0/q0 (%)=', (np.abs(sigqzrec[0]/qzrec[0]))*100.
         
        #======== saving the hz, h'z, qz and om/Lm reconstructions in separated text files 
        savetxt('rec_hz_'+filename+'.dat',rec)
        savetxt('rec_dhz_'+filename+'.dat',drec)
        savetxt('rec_qz_'+filename+'.dat', np.transpose([zrec,qzrec,sigqzrec]))    
        savetxt('rec_omz_lmz_'+filename+'.dat', np.transpose([zrec,om,sigom,Lm/(1.+zrec)**6.,sigLm/(1+zrec)**6.]))   

        #############################################################################################################
        #================ PLOTTING RESULTS ##########################################################################
        #############################################################################################################
        
        ################## plotting the reconstructed hz results
        
        #========= LaTeX rendering text fonts
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

        #=========== Create figure size in inches
        fig, ax1 = plt.subplots(figsize = (10., 7.))
        
        # ========= Define axes
        ax1.set_xlabel(r"$z$", fontsize=22)
        ax1.set_ylabel(r"$H(z)$ ($\mathrm{km} \; \mathrm{s}^{-1} \; \mathrm{Mpc}^{-1}$)", fontsize=22)
        plt.xlim(zmin, zmax+0.02)
        for t in ax1.get_xticklabels(): t.set_fontsize(22)
        for t in ax1.get_yticklabels(): t.set_fontsize(22)
        
        # ========== Plotting the Hz measurements along with the reconstructed H(z) curve with 1, 2 and 3sigma uncertainties 
        plt.errorbar(Z, Hz, yerr=Sigma, fmt='o', color='black')
        ax1.fill_between(zrec, hzrec+1.*sighzrec, hzrec-1.*sighzrec, facecolor='#F08080', alpha=0.80, interpolate=True)
        ax1.fill_between(zrec, hzrec+2.*sighzrec, hzrec-2.*sighzrec, facecolor='#F08080', alpha=0.50, interpolate=True)
        ax1.fill_between(zrec, hzrec+3.*sighzrec, hzrec-3.*sighzrec, facecolor='#F08080', alpha=0.30, interpolate=True)
        plt.legend((r"$H(z)$ rec ($1\sigma$)", "$H(z)$ rec ($2\sigma$)", "$H(z)$ rec ($3\sigma$)", "$H(z)$ data"), fontsize='22', loc='best')
        plt.show()
        
        #========== saving the plot
        fig.savefig('rec_'+filename+'.png')
        
        ################## plotting the reconstructed qz results
        
        #========= LaTeX rendering text fonts
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

        #=========== Create figure size in inches
        fig, ax = plt.subplots(figsize = (10., 7.))

        #========== Define axes
        ax.set_xlabel(r"$z$", fontsize=22)
        ax.set_ylabel(r"$q(z)$", fontsize=22)
        plt.xlim(zmin,zmax+0.02)
        plt.ylim(-1.,1.)
        for t in ax.get_xticklabels(): t.set_fontsize(20)
        for t in ax.get_yticklabels(): t.set_fontsize(20)

        #========= plotting the qz results with 1, 2 and 3sigma uncertainties 
        plt.axhline(0., color='black')
        plt.plot(zrec,qzrec, '--', color='black')
        ax.fill_between(zrec, qzrec+1.*sigqzrec, qzrec-1.*sigqzrec, facecolor='#F08080', alpha=0.80, interpolate=True)
        ax.fill_between(zrec, qzrec+2.*sigqzrec, qzrec-2.*sigqzrec, facecolor='#F08080', alpha=0.50, interpolate=True)
        ax.fill_between(zrec, qzrec+3.*sigqzrec, qzrec-3.*sigqzrec, facecolor='#F08080', alpha=0.30, interpolate=True)
        plt.legend((r"$\mathrm{no \; acc.}$", "$q(z)$", "$1\sigma$", "$2\sigma$", "$3\sigma$"), fontsize='22', loc='best')
        plt.show()
        
        #========== saving the plot
        fig.savefig('rec_qz_'+filename+'.png')
        
        ################## plotting the reconstructed om results
        
        #========= LaTeX rendering text fonts
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

        #=========== Create figure size in inches
        fig, ax2 = plt.subplots(figsize = (10., 7.))
        
        # ========= Define axes
        ax2.set_xlabel(r"$z$", fontsize=22)
        ax2.set_ylabel(r"$\mathcal{O}_\mathrm{m}(z)$", fontsize=22)
        plt.xlim(np.min(Z), zmax+0.02)
        plt.ylim(0.1,0.5)
        for t in ax2.get_xticklabels(): t.set_fontsize(22)
        for t in ax2.get_yticklabels(): t.set_fontsize(22)
        
        # ========== Plotting the om results with 1, 2 and 3sigma uncertainties 
        wm_lcdm = 0.3166
        plt.axhline(wm_lcdm, color='black')
        ax2.fill_between(zrec, om+1.*sigom, om-1.*sigom, facecolor='#F08080', alpha=0.80, interpolate=True)
        ax2.fill_between(zrec, om+2.*sigom, om-2.*sigom, facecolor='#F08080', alpha=0.50, interpolate=True)
        ax2.fill_between(zrec, om+3.*sigom, om-3.*sigom, facecolor='#F08080', alpha=0.30, interpolate=True)
        plt.legend((r"$\Lambda\mathrm{CDM}$", "$\mathcal{O}_\mathrm{m}(z)$ ($1\sigma$)", "$\mathcal{O}_\mathrm{m}(z)$ ($2\sigma$)", "$\mathcal{O}_\mathrm{m}(z)$ ($3\sigma$)"), fontsize='22', loc='best')
        plt.show()
        
        #========== saving the plot
        fig.savefig('rec_om_'+filename+'.png')
        
        ################## plotting the reconstructed Lm results
        
        #========= LaTeX rendering text fonts
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

        #=========== Create figure size in inches
        fig, ax3 = plt.subplots(figsize = (10., 7.))
        
        # ========= Define axes
        ax3.set_xlabel(r"$z$", fontsize=22)
        ax3.set_ylabel(r"$\widehat{\mathcal{L}}_\mathrm{m}(z)$", fontsize=22)
        plt.xlim(np.min(Z), zmax+0.02)
        #plt.ylim(-0.20,0.15)
        for t in ax3.get_xticklabels(): t.set_fontsize(22)
        for t in ax3.get_yticklabels(): t.set_fontsize(22)
        
        # ========== Plotting the om results with 1, 2 and 3sigma uncertainties 
        lm = Lm/(1.+zrec)**6.
        siglm = sigLm/(1.+zrec)**6. 
        plt.axhline(0., color='black')
        ax3.fill_between(zrec, lm+1.*siglm, lm-1.*siglm, facecolor='#F08080', alpha=0.80, interpolate=True)
        ax3.fill_between(zrec, lm+2.*siglm, lm-2.*siglm, facecolor='#F08080', alpha=0.50, interpolate=True)
        ax3.fill_between(zrec, lm+3.*siglm, lm-3.*siglm, facecolor='#F08080', alpha=0.30, interpolate=True)
        plt.legend((r"$\Lambda\mathrm{CDM}$", "$\widehat{\mathcal{L}}_\mathrm{m}(z)$ ($1\sigma$)", "$\widehat{\mathcal{L}}_\mathrm{m}(z)$ ($2\sigma$)", "$\widehat{\mathcal{L}}_\mathrm{m}(z)$ ($3\sigma$)"), fontsize='22', loc='best')
        plt.show()
        
        #========== saving the plot
        fig.savefig('rec_lm_'+filename+'.png')
