#!/usr/bin/env python
# coding: utf-8

# In[1]:


#========
# Author: Maitrey 20220621
# Modification 1: Baibhav 20220626
#========
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import QTable

cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.27, Tcmb0 = 2.725)


# In[6]:


# importing void and galaxy data

data_gal = np.loadtxt('maglim_voidgals_dr7_cbp_102709.dat')
data_void = np.loadtxt('public_void_catalog.dat')

# stripping the data

rdshft_gal = data_gal[:, 2] #redshift of void galaxies
cmvdist_gal = data_gal[:, 3] #comoving distance of void galaxies

x_void = data_void[:, 2] 
y_void = data_void[:, 3]
z_void = data_void[:, 4]
void_radius = data_void[:, 5] 
comoving_dist_void  = np.sqrt(x_void**2+y_void**2+z_void**2)*0.7 #in Mpc h^-1, *h for Mpc

x_gal = data_gal[:, 6]
y_gal = data_gal[:, 7]
z_gal = data_gal[:, 8]
comoving_dist_gal  = np.sqrt(x_gal**2+y_gal**2+z_gal**2)*0.7


# In[10]:


#storing the ra dec
ra_void = data_void[:, 0]
dec_void = data_void[:, 1]
ra_gal = data_gal[:, 0]
dec_gal = data_gal[:, 1]

for i in range(10):#range(len(coord_void)):
    c2 = SkyCoord(ra=ra_gal*u.degree, dec=dec_gal*u.degree, distance=comoving_dist_gal*u.Mpc)
    c1 = SkyCoord(ra=ra_void[i]*u.degree, dec=dec_void[i]*u.degree, distance=comoving_dist_void[i]*u.Mpc)
    dist = c1.separation_3d(c2)
    #dist1 would be of the size of Ngal
    #voidGalxies.append((i+1, j, dist1, void_radius[i]))

    print(np.size(dist)==np.size(ra_gal))


# In[ ]:




