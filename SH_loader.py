#!/usr/bin/env python
import numpy as np
import sys
import time
import os
import matplotlib.pyplot as plt
from pyshtools import SHExpandDH #for expanding a model in spherical harmonic coefficients, on a grid
import sys


class SHloader:
    def __init__(self,datafilepath=None):
        
        #Set up the global grid for model extraction - must be same dimensions as the existing grid text file
        
        self.r_earth = 6371.
        self.nr = 290
        self.nlat = 180
        self.nlon = 360
        self.nl = self.nlat/2 #number of harmonics we're going to end up with
        dlat = 180./self.nlat
        self.depths = np.linspace(0,2890,self.nr)
        self.radii = self.r_earth - self.depths
        
        #Cannot extend to the poles - this way, we make a grid from 0.5 to 89.5 and 0.5 to 359.5 
        self.lons = np.arange(0+dlat/2,360,dlat)
        self.lats = np.arange(-90+dlat/2,90,dlat)
        
        self.datafile = datafilepath
        self.tmpfilename = 'semucb_tmp.npy'
        
        if not os.path.isdir('coeffs'):
            os.system('mkdir coeffs')
        
        if self.datafile:
            print self.datafile
            if os.path.isfile(self.datafile):
                print 'Found input file'
            else:
                print 'Could not find input file: Looking for tmp file'
                
        else:
            if not os.path.isfile(self.tmpfilename):
                print 'Error: no tmp data filename found either: exiting'
                sys.exit(1)
            
    def readtextfile(self):
        
        print 'No tmp file: rereading from model source'
                
        timer = time.time()
        textdata = np.loadtxt(self.datafile)
        np.save(self.tmpfilename,textdata)
        timer = time.time() - timer
        print 'Written tmp model data file: %g s' %timer
    
    def readtmpfile(self):
        
        self.data = np.load(self.tmpfilename)
        print np.shape(self.data)
        self.datashaped = self.data.reshape(self.nr,self.nlat,self.nlon,4)
        print np.shape(self.datashaped)
        #sys.exit(1)
        
    def writecoeffsperdepth(self):
        '''Loop over the radius vector and write out the SH coefficents, one file per depth'''

        data = self.datashaped

        for ir in range(self.nr):
            print 'Writing layer for depth %g' %(int(self.radii[ir]))
            file_lndvs = open('coeffs/lnvs_{:d}.dat'.format(int(self.radii[ir])),'w')
            file_xi = open('coeffs/xi_{:d}.dat'.format(int(self.radii[ir])),'w')

            layer_vs  = np.copy(data[ir,:,:,0]) #Voigt average Vs?
            layer_xi  = np.copy(data[ir,:,:,1]) # Xi is Vsh**2/Nsv**2 ?
            layer_dvs = np.copy(data[ir,:,:,2]) #Relative VA shear velocity
            layer_dxi = np.copy(data[ir,:,:,3]) #Xi perturbation

            layer_lndvs = layer_dvs  #divide by mean for that particular layer

            try:
                np.asarray_chkfinite(layer_lndvs) #Check that all entries are finite - will raise error if not so
                coeffs_lndvs = SHExpandDH(layer_lndvs,sampling=2)
                #sys.exit(1)
                coeffs_xi    = SHExpandDH(layer_xi,sampling=2)
            except ValueError:
                coeffs_lndvs = np.zeros( (2,self.nl,self.nl) )
                coeffs_xi    = np.zeros( (2,self.nl,self.nl) )

            #print np.shape(coeffs_lndvs)

            #---- output lndvs ----
            for l in range(self.nl):
                coeffsa = ['{:4.4e}'.format(c) for c in coeffs_lndvs[0,l,:l+1]]
                coeffsb = ['{:4.4e}'.format(c) for c in coeffs_lndvs[1,l,:l+1]]
                file_lndvs.write(' '.join(coeffsa))
                file_lndvs.write('\n')
                file_lndvs.write(' '.join(coeffsb))
                file_lndvs.write('\n')

            #---- output xi ----
            for l in range(self.nl):
                coeffsa = ['{:4.4e}'.format(c) for c in coeffs_xi[0,l,:l+1]]
                coeffsb = ['{:4.4e}'.format(c) for c in coeffs_xi[1,l,:l+1]]
                file_xi.write(' '.join(coeffsa))
                file_xi.write('\n')
                file_xi.write(' '.join(coeffsb))
                file_xi.write('\n')

            file_lndvs.close()
            file_xi.close()

    def writecoeffsmeandeprange(self,depth1,depth2):
        '''Determine the model mean over a given depth range, then output coeffs file associated with that range'''

        data = self.datashaped

        # Determine if the user's depths are in the radius vector 
        rad1 = self.r_earth - depth1
        rad2 = self.r_earth - depth2

        ind1,ind2 = checkinvector(self.radii,rad1,rad2)

        #file_lndvs = open('coeffs/testlnvs_{:d}.dat'.format(int(self.radii[r])),'w')

        file_meanlndvs = open('coeffs/lnvs_%s_%s_MEAN.dat' %(str(int(depth1)),str(int(depth2))),'w')
        file_meanxi = open('coeffs/xi_%s_%s_MEAN.dat' %(str(int(depth1)),str(int(depth2))),'w')

        count = 0

        meanlayer_vs = np.zeros([self.nlat,self.nlon])
        meanlayer_xi = np.zeros([self.nlat,self.nlon])
        meanlayer_dvs = np.zeros([self.nlat,self.nlon])
        meanlayer_dxi = np.zeros([self.nlat,self.nlon])

        for r in range(ind1,ind2+1):

            layer_vs = np.copy(data[r,:,:,0]) #Voigt average Vs
            layer_xi = np.copy(data[r,:,:,1]) #Xi is Vsh**2/Vsv**
            layer_dvs = np.copy(data[r,:,:,2]) #Relative Voigt average vel 
            layer_dxi = np.copy(data[r,:,:,3]) #Relative Xi perturbation

            meanlayer_dvs = meanlayer_dvs + layer_dvs
            meanlayer_vs = meanlayer_vs + layer_vs
            meanlayer_xi = meanlayer_xi + layer_xi
            meanlayer_dxi = meanlayer_dxi + layer_dxi
            count += 1

        meanlayer_dvs = meanlayer_dvs/count
        meanlayer_vs = meanlayer_vs/count
        meanlayer_xi = meanlayer_xi/count
        meanlayer_dxi = meanlayer_dxi/count

        meanlayer_lndvs = meanlayer_dvs #Not sure why this is here 

        try:
             np.asarray_chkfinite(meanlayer_lndvs) #we don't want to divide by zero

             coeffs_lndvs = SHExpandDH(meanlayer_lndvs,sampling=2)
             coeffs_xi = SHExpandDH(meanlayer_xi,sampling=2)

        except ValueError:

             coeffs_lndvs = np.zeros( (2,self.nl,self.nl) )
             coeffs_xi = np.zeros( (2,self.nl,self.nl) )


        #---- output lndvs ----
        for l in range(self.nl):
            coeffsa = ['{:4.4e}'.format(c) for c in coeffs_lndvs[0,l,:l+1]]
            coeffsb = ['{:4.4e}'.format(c) for c in coeffs_lndvs[1,l,:l+1]]
            file_meanlndvs.write(' '.join(coeffsa))
            file_meanlndvs.write('\n')
            file_meanlndvs.write(' '.join(coeffsb))
            file_meanlndvs.write('\n')

        #---- output xi ----
        for l in range(self.nl):
            coeffsa = ['{:4.4e}'.format(c) for c in coeffs_xi[0,l,:l+1]]
            coeffsb = ['{:4.4e}'.format(c) for c in coeffs_xi[1,l,:l+1]]
            file_meanxi.write(' '.join(coeffsa))
            file_meanxi.write('\n')
            file_meanxi.write(' '.join(coeffsb))
            file_meanxi.write('\n')
        
        file_meanlndvs.close()
        file_meanxi.close()

def checkinvector(invec,d1,d2):
    '''Check if d1 and d2 are in invec: If so, return their indices. If not return the indices of the values closest to them'''

    if d2 > d1:
        print 'Error: d1 > d2!' #because of conversion to radius
        sys.exit(1)

    if d1 in invec:
        occurences = np.where(invec==d1)
        indexd1 = occurences[0][0]
    else:
        invectemp = abs(invec-d1)
        occurences = np.where(invectemp==min(invectemp))
        indexd1 = occurences[0][0]
    if d2 in invec:
        occurences = np.where(invec==d2)
        indexd2 = occurences[0][0]
    else:
        invectemp = abs(invec-d2)
        occurences = np.where(invectemp==min(invectemp))
        indexd2 = occurences[0][0] 

    return indexd1,indexd2     

         
         
if __name__ == '__main__':
   #Testing

   test = SHloader('/Users/rmartinshort/Documents/Berkeley/ASPECT/Tomo_filtering/Harmonic_analysis/semucb_grid.dat.gz')

   #Read the text file (one time use)
   #test.readtextfile()

   #Read the temp file with numpy
   test.readtmpfile()

   #Write the coefficients, one file per depth level
   test.writecoeffsperdepth()

   #Write coefficents for mean levels
   test.writecoeffsmeandeprange(90,140)
   test.writecoeffsmeandeprange(100,200)
   test.writecoeffsmeandeprange(100,400)
   test.writecoeffsmeandeprange(100,600)







