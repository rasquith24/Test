import numpy as np

h = 0.677

#%%

# Defines edges of mass or redshift bins
def make_bin(Nmin, Nmax, dN):
    
    # Uses min, max and bin size to calculate number of bins 
    Nnum = ((Nmax - Nmin) / dN) + 1
    
    # Uses linspace to create equally spaced bins
    Nbin = np.linspace(Nmin, Nmax, num=Nnum)
    
    # Finds the centre of the bins
    Nbinmid = ((Nbin[:-1] + Nbin[1:]) / 2)
    
    # Returns the bin edges and number of bins
    return(Nbin, Nbinmid, Nnum - 1)
    

# Returns the number of galaxies in each bin for the whole sample, passive and star-forming
def N_all(M, z, SFR, Mbin, zbin):
    
    # Corrects for factor of h in SAM output
    M = M / h
    
    # Calculates SSFR and finds galaxies below cut
    SSFR = SFR / M
    SSFRcut = 10**(-11)
    Ip = (SSFR < SSFRcut)
    Mp = M[Ip]
    zp = z[Ip]
    
    # Bins whole sample and passive sample. Finds number in each bin in star-forming sample
    N, Mbin, zbin = np.histogram2d(np.log10(M), z, bins=(Mbin, zbin))
    Np, Mbin, zbin = np.histogram2d(np.log10(Mp), zp, bins=(Mbin, zbin))
    Nsf = N - Np
    
    return(N, Np, Nsf)
    

# Creates the comoving volume array by checking how many redshift slices are in each bin
def comov_vol(z, zbin):
    
    # Finds unique z values which correspond to redshift of slices returned
    zuni = np.unique(z)
    
    # Bins unique values in zbins to see how many slices are in each bin
    Nz, zbin = np.histogram(zuni, bins=zbin)
    
    # Size of the box in comoving Mpc
    box = 125
    
    # Comoving volume in each redshift bin
    CV = Nz * (box / h)**3
    
    return(CV)
    

# Returns the mass function given the number in each bin, the comoving volume and mass bin size
def phi_func(N, CV, dM):
    
    # Mass function - number of galaxies per redshift bin per comoving volume
    phi = np.log10(N / (CV * dM))
    
    return(phi)
    

# Returns the error on the mass function
def phi_err_func(N):
    
    # Relative log error on phi using root N
    phi_rel_err = 0.434 * (1 / np.sqrt(N))
    
    return(phi_rel_err)
    

# Creates an array of masses for each redshift slice depending on where the Schechter fit is trusted to
def Sch_mass(Mmax, Mnum, znum):
    
    # Lower limit of Schechter function fit
    Mlim = np.array([8.37,8.92,9.48,10.03,10.54,10.76])
    
    M = np.zeros((Mnum, znum))
    
    for i in range(0,znum):
        
        # Creates a mass array for each redshift bin from Mlim to Mmax
        M[:,i] = np.linspace(Mlim[i], Mmax, num=Mnum)
        
    return(M)
    

# Returns the value of a single Schechter function
def sin_Sch(M, Ms, phi1, a1):
    
    # Equation for a single Schechter function
    phi = np.log10(phi1 * np.log(10) * ((10**(M - Ms))**(1 + a1)) * np.exp(-10**(M - Ms)))
    
    return(phi)
    

# Returns the value of a double Schechter function
def dou_Sch(M, Ms, phi1, a1, phi2, a2):
    
    # Equation for a double Schechter function
    phi = np.log10(np.log(10) * np.exp(-10**(M - Ms)) * (10**(M - Ms)) * ((phi1 * 10**((M - Ms) * a1)) + (phi2 * 10**((M - Ms) * a2))))
    
    return(phi)
    

# Returns the value of phi for the whole sample
def phi_all_Sch(M, Mnum, znum):
    
    # Values of M*, phi1*, alpha1, phi2* and alpha2 for whole sample
    Ms = np.array([10.97, 11.00, 10.87, 10.81, 10.81, 11.03])
    phi1 = np.array([16.27, 16.25, 13.91, 10.13, 4.79, 1.93]) * 10**(-4)
    a1 = np.array([-0.53, -1.17, -1.02, -0.86, -0.55, -1.01])
    phi2 = np.array([9.47, 0, 0, 0, 0, 0]) * 10**(-4)
    a2 = np.array([-1.37, 0, 0, 0, 0, 0])
    
    phi = np.zeros((Mnum, znum))
    
    for i in range(0, znum):
        
        if i == 0:
            
            # In lowest redshift bin use a double Schechter function
            phi[:,i] = dou_Sch(M[:,i], Ms[i], phi1[i], a1[i], phi2[i], a2[i])
        
        else:
            
            # In other redshift bins use a single Schechter function
            phi[:,i] = sin_Sch(M[:,i], Ms[i], phi1[i] ,a1[i])
            
    return(phi)
    

# Returns the value of phi for the passive sample
def phi_p_Sch(M, Mnum, znum):
    
    # Values of M*, phi1*, alpha1, phi2* and alpha2 for passive sample
    Ms = np.array([10.92, 10.84, 10.73, 10.67, 10.87, 10.80])
    phi1 = np.array([19.68, 14.55, 8.81, 4.15, 1.02, 0.65]) * 10**(-4)
    a1 = np.array([-0.38, -0.36, -0.17, 0.03, -0.71, -0.39])
    phi2 = np.array([0.58, 0.005, 0, 0, 0, 0]) * 10**(-4)
    a2 = np.array([-1.52, -2.32, 0, 0, 0, 0])

    phi = np.zeros((Mnum, znum))
    
    for i in range(0, znum):
        
        if 0 <= i <= 1:
            
            # In two lowest redshift bin use a double Schechter function
            phi[:,i] = dou_Sch(M[:,i], Ms[i], phi1[i], a1[i], phi2[i], a2[i])
        
        else:
            
            # In other redshift bins use a single Schechter function
            phi[:,i] = sin_Sch(M[:,i], Ms[i], phi1[i], a1[i])
            
    return(phi)
    

# Returns the value of phi for the star-forming sample
def phi_sf_Sch(M, Mnum, znum):
    
    # Values of M*, phi* and alpha for star-forming sample
    Ms = np.array([10.81, 10.78, 10.76, 10.85, 10.80, 11.06])
    phi1 = np.array([11.35, 12.71, 8.87, 5.68, 3.72, 1.39]) * 10**(-4)
    a1 = np.array([-1.34, -1.26, -1.21, -1.16, -0.53, -1.03])

    phi = np.zeros((Mnum, znum))
            
    for i in range(0, znum):
        
        # Use a single Schechter function for all redshift bins
        phi[:,i] = sin_Sch(M[:,i], Ms[i], phi1[i], a1[i])
        
    return(phi)
    
