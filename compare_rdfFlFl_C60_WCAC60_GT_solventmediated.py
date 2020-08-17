# this python file is used to compare the full, ref and LMF system of rdfFlFl

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
from matplotlib import rc, font_manager
import scipy.integrate as integrate

rpmf_WCAC60_water=np.loadtxt("2WCAC60_water/profile.xvg", unpack=True)
rpmf_WCAC60_GT=np.loadtxt("2WCAC60_GT_new/profile.xvg", unpack=True)
rpmf_WCAC60_GT095_1=np.loadtxt("2WCAC60_GT095/profile_1.xvg", unpack=True)
rpmf_WCAC60_GT095_2=np.loadtxt("2WCAC60_GT095/profile_2.xvg", unpack=True)


pmf_WCAC60_full=interp.InterpolatedUnivariateSpline(rpmf_WCAC60_water[0,:], rpmf_WCAC60_water[1,:])
pmf_WCAC60_ref=interp.InterpolatedUnivariateSpline(rpmf_WCAC60_GT[0,:], rpmf_WCAC60_GT[1,:])
pmf_WCAC60_095_1=interp.InterpolatedUnivariateSpline(rpmf_WCAC60_GT095_1[0,:], rpmf_WCAC60_GT095_1[1,:])
pmf_WCAC60_095_2=interp.InterpolatedUnivariateSpline(rpmf_WCAC60_GT095_2[0,:], rpmf_WCAC60_GT095_2[1,:])

rpmf_C60_water=np.loadtxt("2C60_water/pmf_2C60_water.xvg", unpack=True)
rpmf_C60_GT=np.loadtxt("2C60_GT/pmf_2C60_GT.xvg", unpack=True)

pmf_C60_full=interp.InterpolatedUnivariateSpline(rpmf_C60_water[0,:], rpmf_C60_water[1,:])
pmf_C60_ref=interp.InterpolatedUnivariateSpline(rpmf_C60_GT[0,:], rpmf_C60_GT[1,:])

ru0=np.loadtxt("2C60_GT/umbrella/table_Fl_Fl.xvg",unpack=True)
u0=interp.InterpolatedUnivariateSpline(ru0[0,:], ru0[5,:])

plt.rc('text', usetex=True)
plt.rc('font', **{'family':'serif', 'serif':['Computer Modern Roman']})

kT=300.0/120.0

av_WCAC60_full=integrate.quad(lambda x:pmf_WCAC60_full(x),1.675,1.725)[0]/0.05
av_WCAC60_ref=integrate.quad(lambda x:pmf_WCAC60_ref(x),1.675,1.725)[0]/0.05
av_WCAC60_095_1=integrate.quad(lambda x:pmf_WCAC60_095_1(x),1.675,1.725)[0]/0.05
av_WCAC60_095_2=integrate.quad(lambda x:pmf_WCAC60_095_2(x),1.675,1.725)[0]/0.05

av_C60_full=integrate.quad(lambda x:pmf_C60_full(x),1.675,1.725)[0]/0.05
av_C60_ref=integrate.quad(lambda x:pmf_C60_ref(x),1.675,1.725)[0]/0.05

r=np.linspace(0.95,1.85,100)

rbin=np.arange(0.95, 1.9, 0.05)


rbin[0]=0.96
#rbin=np.insert(rbin,0,0.96)
#rbin=np.insert(rbin,0,0.952)
plt.figure()
ax1 = plt.subplot(111)
ax1.plot(r,pmf_C60_full(r)-av_C60_full-u0(r)/kT,color="black",lw=3,label="C60; full system")
ax1.plot(rbin,pmf_C60_ref(rbin)-av_C60_ref-u0(rbin)/kT,'s',mec='blue', mfc='none', markersize=18, markeredgewidth=2, alpha=0.75,markevery=1, label="C60; GT system")
ax1.plot(rbin,pmf_WCAC60_full(rbin)-av_WCAC60_full-u0(rbin)/kT,color="purple",ls='--',lw=3,dashes=(8, 8),label="WCAC60; full system")
ax1.plot(rbin,pmf_WCAC60_ref(rbin)-av_WCAC60_ref-u0(rbin)/kT,'o',mec='red', mfc='none', markersize=18, markeredgewidth=2, alpha=0.75,markevery=1, label="WCAC60; GT system")
ax1.plot(rbin,(pmf_WCAC60_095_1(rbin)-av_WCAC60_095_1 + pmf_WCAC60_095_2(rbin)-av_WCAC60_095_2)/2 -u0(rbin)/kT,color="magenta",ls='.-',lw=3,dashes=(4, 4), label="WCAC60; GT2 system")
#ax1.plot(r,pmf_ref(r)-av_ref+u1_LMF(r)/kT,color='red',ls='--',lw=3,dashes=(8, 8),label="SSM")
#ax1.plot(rbin,pmf_ref(rbin)-av_ref+dw(rbin)/kT-dw_av/kT,'o',mec='magenta', mfc='none', markersize=18, markeredgewidth=2, alpha=0.75, label="SSM-no linear")

legend=ax1.legend(bbox_to_anchor=(1.02, 0.535),prop={'size':19.})
frame=legend.get_frame()
ax1.tick_params(axis='both', which='major',length=5, pad=12,top='off', right='off')
plt.xticks( [1.0,1.2,1.4,1.6],fontsize = 30)
plt.yticks( [-12,-8,-4,0,4],fontsize = 30)
plt.xlabel(r"$r$(nm)",fontsize=50)
plt.ylabel(r"$\beta \bar{\omega}_{\mathrm{FlFl}}(r)$", fontsize=50)
plt.xlim(0.935,1.6)
ax1.set_aspect(0.03)
#plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().tight_layout()
plt.savefig("pmf_Fl_WCAFl_full_GT_solventmediated.pdf")


plt.show()