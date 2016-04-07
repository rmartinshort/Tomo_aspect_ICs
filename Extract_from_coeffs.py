#!/usr/bin/env python 

#Use coeff files created by SHLoader to extract grids in S velocity peturbation per depth, and write each grid to a file, which should be exactly the same
#in format at the out.dat file. This will then be used by Extraction scripts. 

import sys
import os
from pyshtools import MakeGridDH
import glob
import numpy as np
import matplotlib.pyplot as plt


def ArrangeCoeffFiles(coeffsdir):
	'''Enter coeffs directory and make a sorted list of the depths''' 

	if os.path.exists(coeffsdir):

		coeffiles = []

		os.chdir(coeffsdir)
		lnvs = glob.glob('lnvs_*.dat')
		for coeffile in sorted(lnvs):

			if 'MEAN' not in coeffile:
				coeffiles.append(coeffile)

		return coeffiles

def read_coeffs(fname,nl):

    infile = open(fname,'r')
    coeffs = np.zeros( (2,nl,nl) )

    for l in range(nl):
        alm = infile.readline().strip().split()
        blm = infile.readline().strip().split()
        coeffs[0,l,:l+1] = [float(c) for c in alm]
        coeffs[1,l,:l+1] = [float(c) for c in blm]
    return coeffs


def Taper_harmonic_coeffs(coeffs,llmin,llmax):
	'''Apply a simple cosine taper to the spherical harmonics list'''

	#Lets have the taper span 20 degrees, so we start from lmax-6 and end at lmax+6

	#####First taper at the 'maximum' end 
	thw = 20 #half width of the cosine taper, in degrees

	if llmax < 80:

		print 'Tapering llmax with half width: %g' %thw

		#any coeffs beyond the upper limit of the taper are set to zero
		coeffs[:,llmax+thw+1:,:] = 0
		taperstart = llmax - thw

		#apply the taper to all coeffs between llmax-thw and llmax+thw
		for coeffnumber in range(0,thw*2):
			cf = taperstart + coeffnumber

			#print 'Coeff_number, multiplying factor'
			#print cf, np.cos(2*np.pi*coeffnumber/(thw*8))

			A0 = coeffs[:,cf,:].copy()
			coeffs[:,cf,:] = A0*np.cos(2*np.pi*coeffnumber/(thw*8))

	else:
		print 'Trying to taper to coefficients beyond the length of the array!'
		sys.exit(1)


	#####Now taper at the minumum end, if necassary
	if llmin > 30:

		print 'tapering llmin with half width: %g' %thw

		#any coeffs beyond the upper limit of the taper are set to zero
		coeffs[:,:llmin-thw,:] = 0
		taperstart = llmin + thw
		print llmin
		print taperstart

		#apply the taper to all coeffs between llmax-thw and llmax+thw
		for coeffnumber in range(0,thw*2):
			cf = taperstart - coeffnumber

			#print 'Coeff_number, multiplying factor'
			#print cf, np.cos(2*np.pi*coeffnumber/(thw*8))

			A0 = coeffs[:,cf,:].copy()
			coeffs[:,cf,:] = A0*np.cos(2*np.pi*coeffnumber/(thw*8))

	else:
		print 'Trying to taper to coefficients beyond the bottom of the array!\n Skipping bottom taper'


	return coeffs

def Append2master(openmasterfile,coeffsdir,lmin,lmax):
	'''Loop though ordered coefficient files, cut out coefficients outside the provided range, make global grid and append points to the 'master' file'''

	coeffiles = ArrangeCoeffFiles(coeffsdir)

	outfile = openmasterfile

    #Want all but the last 7 coeffs files, which corresponds to the model at the surface so is probably incorrect ( d < 70km )
	for coeffile in coeffiles[:-7]:

		path = coeffsdir+'/'+coeffile
		print 'Extracting from %s' %path

		rad = float(coeffile.split('_')[1].split('.')[0])
		print 'Dealing with radius: %g' %rad

		coeffsuntapered = read_coeffs(coeffile,nl=90)

		#Window just the selected coeffs (apply cosine taper)
		#coeffs[:,:lmin,:] = 0
		#coeffs[:,lmax+1:,:]=0
		coeffs = Taper_harmonic_coeffs(coeffsuntapered,lmin,lmax)

		#make global grid
		grid  = MakeGridDH(coeffs,sampling=2)

		#Append global grid points to the masterfile
		Gshape = np.shape(grid)
		Glat = Gshape[0]
		Glon = Gshape[1]

		grid = np.transpose(grid)

		#plt.imshow(np.transpose(grid))
		#plt.xlabel('x')
		#plt.ylabel('y')
		#plt.show()
		#sys.exit(1)

		#loop over lat and lon grids: lon has to vary before lat

		for i in range(Glat):
		 	for j in range(Glon):

		 		lon =  j
		 		lat = -90 + i
		 		val = grid[j,i]
		 		outfile.write('%g %g %g %g\n' %(rad,lon,lat,val))

		 	outfile.write('%g 360 %g %g\n' %(rad,lat,val))

		#Add extra line for the N pole
		for i in range(Glon+1):
			outfile.write('%g %g 90 %g\n' %(rad,i,val))

	outfile.close()



if __name__ == '__main__':

	#ArrangeCoeffFiles('../coeffs')

	#Open a file yo output the results to
	out = open('test_L1_L20.dat','w')

	#Run the main functtion to generate the grids from these coefficients (from llmax to llmin)
	Append2master(out,'../coeffs',0,20)



