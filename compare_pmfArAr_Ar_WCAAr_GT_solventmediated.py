# this python file is used to compare the full, ref and LMF system of rdfArAr

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
from matplotlib import rc, font_manager
import scipy.integrate as integrate

rpmf_Ar_water=np.loadtxt("2Ar_water/pmf_2Ar_water.xvg", unpack=True)
rpmf_Ar_GT=np.loadtxt("2Ar_GT/pmf_2Ar_GT.xvg", unpack=True)

pmf_Ar_full=interp.InterpolatedUnivariateSpline(rpmf_Ar_water[0,:], rpmf_Ar_water[1,:])
pmf_Ar_ref=interp.InterpolatedUnivariateSpline(rpmf_Ar_GT[0,:], rpmf_Ar_GT[1,:])

rpmf_WCAAr_water=np.loadtxt("2WCAAr_water/profile.xvg", unpack=True)
rpmf_WCAAr_GT=np.loadtxt("2WCAAr_GT/profile.xvg", unpack=True)

pmf_WCAAr_full=interp.InterpolatedUnivariateSpline(rpmf_WCAAr_water[0,:], rpmf_WCAAr_water[1,:])
pmf_WCAAr_ref=interp.InterpolatedUnivariateSpline(rpmf_WCAAr_GT[0,:], rpmf_WCAAr_GT[1,:])

ru0=np.loadtxt("table_Ar_Ar.xvg",unpack=True)
u0=interp.InterpolatedUnivariateSpline(ru0[0,:], ru0[5,:])

plt.rc('text', usetex=True)
plt.rc('font', **{'family':'serif', 'serif':['Computer Modern Roman']})

kT=300.0/120.0

av_Ar_full=integrate.quad(lambda x:pmf_Ar_full(x),1.175,1.225)[0]/0.05
av_Ar_ref=integrate.quad(lambda x:pmf_Ar_ref(x),1.175,1.225)[0]/0.05

av_WCAAr_full=integrate.quad(lambda x:pmf_WCAAr_full(x),1.175,1.225)[0]/0.05
av_WCAAr_ref=integrate.quad(lambda x:pmf_WCAAr_ref(x),1.175,1.225)[0]/0.05

r=np.linspace(0.31,1.2,100)
rbin=np.linspace(0.32,1.2,30)
#rbin=np.insert(rbin,0,0.305)
plt.figure()
ax1 = plt.subplot(111)
ax1.plot(r,pmf_Ar_full(r)-av_Ar_full-u0(r)/kT,color="black",lw=3,label="Ar; full water")
ax1.plot(rbin,pmf_Ar_ref(rbin)-av_Ar_ref-u0(rbin)/kT,'s',mec='blue', mfc='none', markersize=18, markeredgewidth=2, alpha=0.75,markevery=1,label="Ar; GT water")
ax1.plot(r,pmf_WCAAr_full(r)-av_WCAAr_full-u0(r)/kT,color="purple",ls='--',lw=3,dashes=(8, 8),label="WCAAr; full water")
ax1.plot(rbin,pmf_WCAAr_ref(rbin)-av_WCAAr_ref-u0(rbin)/kT,'o',mec='red', mfc='none', markersize=18, markeredgewidth=2, alpha=0.75,markevery=1,label="WCAAr; GT water")
#ax1.plot(r,pmf_ref(r)-av_ref+u1_LMF(r)/kT,color='red',ls='--',lw=3,dashes=(8, 8),label="SSM")
#ax1.plot(rbin,pmf_ref(rbin)-av_ref+dw(rbin)/kT-dw_av/kT+dw_coloumb(rbin)/kT-dw_coloumb_av/kT,'o',mec='magenta', mfc='none', markersize=18, markeredgewidth=2, alpha=0.75, label="SSM-no linear")

legend=ax1.legend(bbox_to_anchor=(1.02, 0.62),prop={'size':25})
frame=legend.get_frame()
ax1.tick_params(axis='both', which='major',length=5, pad=12,top='off', right='off')
plt.xticks( [0.4,0.6,0.8],fontsize = 35)
plt.yticks( [-3, -2, -1, 0, 1], fontsize = 35)
plt.xlabel(r"$r$(nm)",fontsize=50)
plt.ylabel(r"$\beta \bar{\omega}_{\mathrm{ArAr}}(r)$", fontsize=50)
plt.xlim(0.285,0.9)
#plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().tight_layout()
plt.savefig("pmf_Ar_WCAAr_GT_solventmediated.pdf")


plt.show()