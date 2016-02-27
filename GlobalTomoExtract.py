#!/usr/bin/env python

#RMS Feb 2016
#Tidied version of ExtractGlobalTomo.py, which can extract from a global tomographic database to ASPECT input
#Will extract absolute velocities or perturbations, athough right now perturbations are preferrable

import argparse
import numpy as np
from evtk.hl import pointsToVTK
from scipy.interpolate import RegularGridInterpolator
import os
import sys

#Import functions to generate absolute temperature perturbation
import Absolute_Temps as AT


def main():

	'''

	Workflow is as follows
	> Read data from tomo input file on format rad lon lat V1 V2 etc, where lon varies before lat and rad varies last
	> Use this data to make a 3D interpolation function
	> Construct a grid onwhich to output the values in a format ASPECT can read
	> Interpolate the tomography values on this grid, and covnvert to absolute temperature using the functions in Absolute_Temps
	> Sweep the output file and convert to colatitude
	> Write a .vtu file

	'''

	parser = argparse.ArgumentParser()

	parser.add_argument('-extract',action='store',dest='extracttype',help='Choose the variable to extract: Choose from SH SV SM PertS or HPertS, for PertS components from the harmonic expansiod')
	parser.add_argument('-datafile',action='store',dest='userdatafile',help='Datafile from which to extract tomography, if present. If not, append NONE and one will be created using Scotts extraction code')

	userinput = parser.parse_args()

	#Basename for the output file
	OutfileID='ASPECT_tomoglobe_absolutetemps_L1_L5'

	if userinput.userdatafile == 'NONE':
		extractfromdataset()
		datafile = 'out.dat'
	else:
		if not os.path.exists(userinput.userdatafile):
			print 'Chosen datafile %s does not appear to exist in this directory' %datafile

	datafile=open(userinput.userdatafile,'r')
	lines = datafile.readlines()
	datafile.close()

    #Value to multiply the output grid by (only not equal to 1 when we extract from a file made by the spherical harmonic code, which 
    #makes values as fractions rather than percentages)
	M = 1 

	if userinput.extracttype == 'PertS':
		fcount = 3
	elif userinput.extracttype == 'SM':
		fcount = 'Mean'
	elif userinput.extracttype == 'SV':
		fcount = -2
	elif userinput.extracttype == 'SH':
		fcount = -1
	elif userinput.extracttype == 'HPertS':
		fcount = -1
		M = 100
	else:
		print 'User input data type not recognized!'
		sys.exit(1)


	#-----------------------------------------------------------
	#Generate interpolation function using the input grid 
	#-----------------------------------------------------------
	interpfunc,LATGRID = Geninterpgrid(lines,fcount,Mulval=M)

	#-----------------------------------------------------------
	#Now make out grids and interpolate using the provided function 
	#-----------------------------------------------------------

	#These values the output grid: We can have the same number of points at the original grid, but this is not necessary

	minrad = 3490.0

	#We can't extract datapoints above 100km from the surface, as these are probably anomalous: However, we still need a rule for this layer, so we have to write the
	#if statements seen in the triple loop below

	maxrad = 6370.0 

	latstepnumber = 181
	lonstepnumber = 361
	radstepnumber = 121

	radgrid = np.linspace(minrad,maxrad,radstepnumber)*1000
	latgrid = np.linspace(-90,90,latstepnumber)
	longrid = np.linspace(0,360,lonstepnumber)

	#----------------------------------------------
	#Interpolate the dataset at the desired points (in the order expected by ASPECT)
	#----------------------------------------------

    #This is a temporary output file
	outfilename = 'ASPECT_tomoglobe_%s.dat' %userinput.extracttype
	outfile = open(outfilename ,'w')

	#generate the meterial property interpolation functions, from Ian's curves 
	f_temps,f_rhos,f_Ks,f_As,f_vs,fmeanvs,f_etas = AT.generateinterpfuncs()

	for lat in latgrid:
	  	for lon in longrid:
	  		for rad in radgrid:

	  			### Here we've assumed that absolute velocity has been extracted: Find the material properties that correspond to this depth value
		  		tr = f_temps(rad)
		  		#trho = f_rhos(rad)
		  		tK = f_Ks(rad)
		  		tA = f_As(rad)
		  		#tv = f_vs(rad)

	  			if rad < 6270000: #Must be less than 100 km from the surface the extract from the interpolation function

					interpoint = interpfunc([rad,lon,lat]) #This should be some measure of velocity

					etaval = f_etas(rad) #get the dln(rho)/dln(V) value for this depth

					#meanvs = fmeanvs(rad) #The global model mean velocity at this depth
					
					#print 'Depth: %g' %(maxrad*1000-rad)
					#print 'model vel: %g' %tv
					
					#print 'interp vel: %g' %interpoint
					#print '-------------------'

					#Method 1: Simple: Don't use shear modulus to determine Rhodiff (seems to be more accurate!)

					#interpoint = 0 #Used if we just want the 3D radial temperature perturbation 
					temp = AT.ReturnTempPert(rad, interpoint, tr, tA,eta=etaval)

					#Method 2: Use the absolute velocoty difference between the actual velocity at some depth and the mean velocity to determin rhodiff 
					#temp = AT.ReturnTempAbsolute(rad, interpoint, meanvs, tr, tK, tA)
				else:
					Svel = 0.0 #estimate of S velocity for the 'crustal' part of the model - really the lithopshere since its anything above 100km. Set at 0.0 so that we get the radial
					#profile geotherm here
					temp = AT.ReturnTempPert(rad, Svel, tr, tA)

	  			outfile.write('%g %g %g %g\n' %(rad,lon,lat,temp))

	outfile.close()

	#-----------------------------------------------------------
	#Now, reopen and reorder so that latitude is converted to colatitude and the order of the points remains the same
	#-----------------------------------------------------------

	datafile=open(outfilename,'r')
	outfilename2  = OutfileID+'_%s.dat' %userinput.extracttype
	outfile = open(outfilename2,'w')

	outfile.write('# POINTS: %g %g %g\n' %(radstepnumber,lonstepnumber,latstepnumber))
	lines = datafile.readlines()
	datafile.close()

	tmpX = []
	tmpY = []
	tmpZ = []
	tmpData = []

	i = 0
	for lat in sorted(90-LATGRD):

		for line in lines:
			vals = line.split()
			radf = float(vals[0])
			latf = 90-float(vals[2])
			lonf = float(vals[1])
			tomof = float(vals[3])
			if (lat==latf):

				tmpX.append(lonf*(np.pi/180.0))
				tmpY.append(latf*(np.pi/180.0))
				tmpZ.append(radf)
				tmpData.append(tomof)

	            #Previous, crude estimation of the temperature 
				#Temporary for december test: Convert to temperature in a crude fashion. Need to link this to a geothermal profile at some point
				#tomotemp = (1600 - ((tomof/100.0)*0.25)/1e-4)

				#if tomotemp < 1300:
				#	tomotemp = 1300
				#if tomotemp > 1800:
				#	tomotemp = 1700

				outfile.write('%g %g %g %g\n' %(radf,lonf*(np.pi/180),latf*(np.pi/180),tomof))
			
		print 'Done sweep %g' %i
		i+=1


	#Make a point cloud file containing the full grid - this allows us to check that the data has been correctly entered 

	X,Y,Z = spheretocart(np.array(tmpX),np.array(tmpY),np.array(tmpZ))
	pointsToVTK("./ASPECT_GlobeTomo_input",X,Y,Z,data = {"Pertb": np.array(tmpData)})



def Geninterpgrid(lines,fcount,Mulval=1):

	'''
	-----------------------------------------------------------
	GENERIC METHOD FOR OBTAINING THE GRID SPACINGS FROM A SPHERICAL GRID FILE
	-----------------------------------------------------------
	Determine the grids in lat, lon and depth, and their spacings 

	'''

	#First, load the data

	DATAVEC = np.zeros(len(lines))
	i=0
	for line in lines:
		vals = line.split()

		if fcount != 'Mean': #This is what we'll use to extract the absolute temperatures
			VAL = float(vals[fcount])
		else:
			VAL = (float(vals[-2]) + float(vals[-1]))/2.0

		DATAVEC[i] = VAL  #*1000 #Convert to m/s
		i+= 1

	#Get lon grid
	tmplon = []
	for line in lines:
		vals = line.split()
		lon = float(vals[1])

		if lon in tmplon:
			break 
		else:
			tmplon.append(lon)

	#Get latitude grid
	tmplat = []
	for line in lines[0:len(lines):len(tmplon)]:
		vals = line.split()
		lat = (float(vals[2]))

		if lat in tmplat:
			break 
		else:
			tmplat.append(lat)

	#Get radius grid
	tmpradius = []
	for line in lines[0:len(lines):len(tmplon)*len(tmplat)]:
		vals = line.split()
		rad = float(vals[0])

		if rad in tmpradius:
			break 
		else:
			tmpradius.append(rad)

	#Determine the number of points in lat,lon and radius
	Rpoints = len(tmpradius)
	Latpoints = len(tmplat)
	Lonpoints = len(tmplon)

	Radiusbounds = [min(tmpradius),max(tmpradius)]
	Latbounds = [min(tmplat),max(tmplat)]
	Lonbounds = [min(tmplon),max(tmplon)]

	#----------------------------------------------
	#Create grid for interpolation
	#----------------------------------------------
	#Generate the grids

	RGRD = np.array(tmpradius)*1000.0 #convert to meters - all input grids are in km 
	LATGRD = np.array(tmplat)
	LONGRD = np.array(tmplon)

	#Data must be a 3D array
	DATA = np.zeros((Rpoints,Lonpoints,Latpoints))

	for i in range(Rpoints):
		for j in range(Latpoints):
			for k in range(Lonpoints):
				counter = i*Latpoints*Lonpoints + j*Lonpoints + k
				DATA[i,k,j] = DATAVEC[counter]*Mulval

	print "==================================================="
	print 'Interpolating on the following grids'
	print "==================================================="
	print 'Latitude:'
	print "==================================================="
	print LATGRD
	print "==================================================="
	print 'Radius'
	print "==================================================="
	print RGRD
	print "==================================================="
	print 'Longitude'
	print "==================================================="
	print LONGRD

	interpfunc = RegularGridInterpolator((RGRD,LONGRD,LATGRD),DATA,method='nearest')

	#LATGRD is needed later in the code

	return interpfunc,LATGRD


def spheretocart(theta,phi,rad):
	#Convert spherical to cartesian coodinates, for creating point-cloud vtu files
	#Input should be LON,COLAT,RAD

	x = rad*np.sin(phi)*np.cos(theta)
	y = rad*np.sin(phi)*np.sin(theta)
	z = rad*np.cos(phi)

	return x,y,z


def extractfromdataset():
	'''Calls Scott's code to extract data from the model at a grid of points, generated here'''

	minrad = 3481.0
	maxrad = 6271.0 #We ignore the upper 100 km, because the tomographic models are not well resolved there. In any case, we're only interested in the mantle flow

	#Set up the extraction grid
	radstepnumber = 121
	latstepnumber = 181
	lonstepnumber = 361


	radgrid = np.linspace(minrad,maxrad,radstepnumber)
	latgrid = np.linspace(-90,90,latstepnumber)
	longrid = np.linspace(0,360,lonstepnumber)

	lonlatfile = open('lonlattest.dat','w')

	for lat in latgrid:
		for lon in longrid:
			lonlatfile.write('%g %g\n' %(lon,lat))
	lonlatfile.close()

	radiusstring = '-r '
	for rad in radgrid:
		radiusstring += '%s -r ' %str(rad)

	radiusstring = radiusstring[:-3]

	if os.path.isfile('./a3d_dist.x'):
		os.system('./a3d_dist.x %s -f lonlattest.dat -o out.dat' %radiusstring)
		os.system('rm lonlattest.dat')

	else:
		print 'Tried to access Scotts code to extract from SEMUBCWM1 database, but cannot find executable a3d_dist.x'


if __name__ == '__main__':

	main()








