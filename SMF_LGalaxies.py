# Plots the Schechter functions from Muzzin et al. (2013) and the mass functions for LGalaxies in 6 redshift bins between 0 < z < 3

#%%
# Importing modules

# Importing numpy and shortening to np
import numpy as np

# Importing matplotlib and shortening to plt
import matplotlib.pyplot as plt

# Importing pyfits to read in fits data
import pyfits

# Importing own function file
from functions import *
    
#%% 
# Reading LGalaxies data in

# Reading LGalaxies fits file in
data = pyfits.getdata('LGalaxies.fits', 1)

# Pulling out columns that are needed
Mass = data.field('M')
SFR = data.field('SFR')
redshift = data.field('z')

#%%
# Defining bin edges and binning data

# Left edge, right edge and bin size for mass and redshift bins
Mmin = 9
Mmax = 12
dM = 0.1
zmin = 0
zmax = 3
dz = 0.5

# Uses function to create bins and passes back edges and number of bins
[Mbin, M, Mnum] = make_bin(Mmin, Mmax, dM)
[zbin, z, znum] = make_bin(zmin, zmax, dz)

# Ensures they are just integers to be used later
Mnum = int(Mnum)
znum = int(znum)

# Bins the data. Also uses a SSFR cut to select passive galaxies 
[N, Np, Nsf] = N_all(Mass, redshift, SFR, Mbin, zbin)

# Calculates the comoving volume in each redshift bin
CV = comov_vol(redshift, zbin)

#%%
# Mass functions and passive fraction for LGalaxies

# Uses the function to calculate the mass function and error for the whole sample, passive sample and star-forming sample
phi = phi_func(N, CV, dM)
phi_err = phi_err_func(N)

phip = phi_func(Np, CV, dM)
phip_err = phi_err_func(Np)

phisf = phi_func(Nsf, CV, dM)
phisf_err = phi_err_func(Nsf)

# Calculates passive fraction and error
pf = Np / N
pf_err = 1 / np.sqrt(N)

#%%
# Schechter functions

# Creates a mass array for the Schechter functions. This changes with redshift as the lower mass limit changes
M_Sch = Sch_mass(Mmax, Mnum, znum)

# Calculates the values of the Schechter functions at M_Sch
phi_Sch = phi_all_Sch(M_Sch, Mnum, znum)
phip_Sch = phi_p_Sch(M_Sch, Mnum, znum)
phisf_Sch = phi_sf_Sch(M_Sch, Mnum, znum)

#%%
# Plotting all redshift bins

# Is there an easier way to plot all these without having to write them all out? Maybe loop through them?

# Number of plots across and down
xnum = 6
ynum = 3

# Defining figure and size
fig = plt.figure(figsize=(3 * xnum, 3 * ynum))

# Adding subplots and setting which axis they share
ax1 = fig.add_subplot(ynum, xnum, 1)
ax2 = fig.add_subplot(ynum, xnum, 2, sharex=ax1, sharey=ax1)
ax3 = fig.add_subplot(ynum, xnum, 3, sharex=ax1, sharey=ax1)
ax4 = fig.add_subplot(ynum, xnum, 4, sharex=ax1, sharey=ax1)
ax5 = fig.add_subplot(ynum, xnum, 5, sharex=ax1, sharey=ax1)
ax6 = fig.add_subplot(ynum, xnum, 6, sharex=ax1, sharey=ax1)
ax7 = fig.add_subplot(ynum, xnum, 7, sharex=ax1, sharey=ax1)
ax8 = fig.add_subplot(ynum, xnum, 8, sharex=ax1, sharey=ax1)
ax9 = fig.add_subplot(ynum, xnum, 9, sharex=ax1, sharey=ax1)
ax10 = fig.add_subplot(ynum, xnum, 10, sharex=ax1, sharey=ax1)
ax11 = fig.add_subplot(ynum, xnum, 11, sharex=ax1, sharey=ax1)
ax12 = fig.add_subplot(ynum, xnum, 12, sharex=ax1, sharey=ax1)
ax13 = fig.add_subplot(ynum, xnum, 13, sharex=ax1)
ax14 = fig.add_subplot(ynum, xnum, 14, sharex=ax1, sharey=ax13)
ax15 = fig.add_subplot(ynum, xnum, 15, sharex=ax1, sharey=ax13)
ax16 = fig.add_subplot(ynum, xnum, 16, sharex=ax1, sharey=ax13)
ax17 = fig.add_subplot(ynum, xnum, 17, sharex=ax1, sharey=ax13)
ax18 = fig.add_subplot(ynum, xnum, 18, sharex=ax1, sharey=ax13)

# Plotting mass function and Schechter function for whole sample for all redshift bins
ax1.errorbar(M, phi[:,0], yerr=phi_err[:,0], fmt='kx')
ax2.errorbar(M, phi[:,1], yerr=phi_err[:,1], fmt='kx')
ax3.errorbar(M, phi[:,2], yerr=phi_err[:,2], fmt='kx')
ax4.errorbar(M, phi[:,3], yerr=phi_err[:,3], fmt='kx')
ax5.errorbar(M, phi[:,4], yerr=phi_err[:,4], fmt='kx')
w = ax6.errorbar(M, phi[:,5], yerr=phi_err[:,5], fmt='kx')

ax1.plot(M_Sch[:,0], phi_Sch[:,0], 'k')
ax2.plot(M_Sch[:,1], phi_Sch[:,1], 'k')
ax3.plot(M_Sch[:,2], phi_Sch[:,2], 'k')
ax4.plot(M_Sch[:,3], phi_Sch[:,3], 'k')
ax5.plot(M_Sch[:,4], phi_Sch[:,4], 'k')
mw, = ax6.plot(M_Sch[:,5], phi_Sch[:,5], 'k')

# Plotting mass function and Schechter function for passive sample for all redshift bins
ax7.errorbar(M, phip[:,0], yerr=phip_err[:,0], fmt='rx')
ax8.errorbar(M, phip[:,1], yerr=phip_err[:,1], fmt='rx')
ax9.errorbar(M, phip[:,2], yerr=phip_err[:,2], fmt='rx')
ax10.errorbar(M, phip[:,3], yerr=phip_err[:,3], fmt='rx')
ax11.errorbar(M, phip[:,4], yerr=phip_err[:,4], fmt='rx')
p = ax12.errorbar(M, phip[:,5], yerr=phip_err[:,5], fmt='rx')

ax7.plot(M_Sch[:,0], phip_Sch[:,0], 'r')
ax8.plot(M_Sch[:,1], phip_Sch[:,1], 'r')
ax9.plot(M_Sch[:,2], phip_Sch[:,2], 'r')
ax10.plot(M_Sch[:,3], phip_Sch[:,3], 'r')
ax11.plot(M_Sch[:,4], phip_Sch[:,4], 'r')
mp, = ax12.plot(M_Sch[:,5], phip_Sch[:,5], 'r')

# Plotting mass function and Schechter function for star-forming sample for all redshift bins
ax7.errorbar(M, phisf[:,0], yerr=phisf_err[:,0], fmt='bx')
ax8.errorbar(M, phisf[:,1], yerr=phisf_err[:,1], fmt='bx')
ax9.errorbar(M, phisf[:,2], yerr=phisf_err[:,2], fmt='bx')
ax10.errorbar(M, phisf[:,3], yerr=phisf_err[:,3], fmt='bx')
ax11.errorbar(M, phisf[:,4], yerr=phisf_err[:,4], fmt='bx')
sf = ax12.errorbar(M, phisf[:,5], yerr=phisf_err[:,5], fmt='bx')

ax7.plot(M_Sch[:,0], phisf_Sch[:,0], 'b')
ax8.plot(M_Sch[:,1], phisf_Sch[:,1], 'b')
ax9.plot(M_Sch[:,2], phisf_Sch[:,2], 'b')
ax10.plot(M_Sch[:,3], phisf_Sch[:,3], 'b')
ax11.plot(M_Sch[:,4], phisf_Sch[:,4], 'b')
msf, = ax12.plot(M_Sch[:,5], phisf_Sch[:,5], 'b')

# Plotting red fraction for all redshift bins
ax13.errorbar(M, pf[:,0], yerr=pf_err[:,0], fmt='kx')
ax14.errorbar(M, pf[:,1], yerr=pf_err[:,1], fmt='kx')
ax15.errorbar(M, pf[:,2], yerr=pf_err[:,2], fmt='kx')
ax16.errorbar(M, pf[:,3], yerr=pf_err[:,3], fmt='kx')
ax17.errorbar(M, pf[:,4], yerr=pf_err[:,4], fmt='kx')
pfp = ax18.errorbar(M, pf[:,5], yerr=pf_err[:,5], fmt='kx')

# Setting the space between plots to 0 and turning axis ticks off
fig.subplots_adjust(hspace=0, wspace=0)
plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)
plt.setp(ax4.get_xticklabels(), visible=False)
plt.setp(ax5.get_xticklabels(), visible=False)
plt.setp(ax6.get_xticklabels(), visible=False)
plt.setp(ax7.get_xticklabels(), visible=False)
plt.setp(ax8.get_xticklabels(), visible=False)
plt.setp(ax9.get_xticklabels(), visible=False)
plt.setp(ax10.get_xticklabels(), visible=False)
plt.setp(ax11.get_xticklabels(), visible=False)
plt.setp(ax12.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)
plt.setp(ax3.get_yticklabels(), visible=False)
plt.setp(ax4.get_yticklabels(), visible=False)
plt.setp(ax5.get_yticklabels(), visible=False)
plt.setp(ax6.get_yticklabels(), visible=False)
plt.setp(ax8.get_yticklabels(), visible=False)
plt.setp(ax9.get_yticklabels(), visible=False)
plt.setp(ax10.get_yticklabels(), visible=False)
plt.setp(ax11.get_yticklabels(), visible=False)
plt.setp(ax12.get_yticklabels(), visible=False)
plt.setp(ax14.get_yticklabels(), visible=False)
plt.setp(ax15.get_yticklabels(), visible=False)
plt.setp(ax16.get_yticklabels(), visible=False)
plt.setp(ax17.get_yticklabels(), visible=False)
plt.setp(ax18.get_yticklabels(), visible=False)

# Setting axis range
ax1.axis([8.8, 12.2, -6.5, -1.5])
ax13.axis([8.8, 12.2, 0, 1])

# Labelling axis and putting on titles
ax13.set_xlabel('log(M [M$_\odot$])')
ax14.set_xlabel('log(M [M$_\odot$])')
ax15.set_xlabel('log(M [M$_\odot$])')
ax16.set_xlabel('log(M [M$_\odot$])')
ax17.set_xlabel('log(M [M$_\odot$])')
ax18.set_xlabel('log(M [M$_\odot$])')
ax1.set_ylabel('log($\phi$ [Mpc$^{-3}$dex$^{-1}$])')
ax7.set_ylabel('log($\phi$ [Mpc$^{-3}$dex$^{-1}$])')
ax13.set_ylabel('Passive fraction')
ax1.set_title('0.0<z<0.5')
ax2.set_title('0.5<z<1.0')
ax3.set_title('1.0<z<1.5')
ax4.set_title('1.5<z<2.0')
ax5.set_title('2.0<z<2.5')
ax6.set_title('2.5<z<3.0')

# Creating subtitle, legend and saving plot
suptitle = plt.suptitle('LGalaxies', fontsize=15)
legent = ['Muzzin et al. (2013) - all', 'Muzzin et al. (2013) - passive', 'Muzzin et al. (2013) - star-forming', 'LGalaxies - all', 'LGalaxies - passive', 'LGalaxies - star-forming', 'LGalaxies - passive fraction']
leg = plt.figlegend([mw, mp, msf, w, p, sf, pfp], legent, loc=(0.81, 0.70), fontsize=10)
fig.subplots_adjust(right=0.8, bottom=0.2)
fig.savefig('LGalaxies_SMF.png', bbox_extra_artists=(leg, suptitle), bbox_inches='tight')

# Showing plot
plt.show()

