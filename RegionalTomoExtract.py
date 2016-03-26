#!/usr/bin/env python

#Convert a regional tomography file (in this case DNA13SVJ) into ASPECT input. This script will load the regional file and generate a global temperature initial condition for ASPECT
#only with pertrubations within the region that is described by the input tomography. At some point it would be useful to load multiple regional tomographies. 

import numpy as np
from math import sqrt
from evtk.hl import pointsToVTK
from scipy.interpolate import RegularGridInterpolator

#Import functions to generate absolute temperature perturbation
import Absolute_Temps as AT

#-----------------------------------------------------------
#Helpful functions
#-----------------------------------------------------------

def converttoASPECT(R,lon,lat):
	'''Converts coordinates in R [km], lon [deg] and lat [deg] into ASPECT required format: R [m], lon [rad], colat [rad]'''

	r = R*1000.0
	rlon = lon*(np.pi/180.0)
	rlat = lat*(np.pi/180.0)

	return r,rlon,rlat


def spheretocart(theta,phi,rad):
	#Convert spherical to cartesian coodinates, for creating point-cloud vtu files
	#Input should be LON,COLAT,RAD

	x = rad*np.sin(phi)*np.cos(theta)
	y = rad*np.sin(phi)*np.sin(theta)
	z = rad*np.cos(phi)

	return x,y,z

def colatconvert(tmplat,infilelines):
	'''Converts an input file to radius longitude colatitude tomoval'''

	outfile = open('exout_colat_tmp.dat','w')

	i = 0
	for lat in sorted(tmplat):

		for line in infilelines:
			vals = line.split('   ')
			radf = float(vals[0])
			latf = 90-float(vals[1])
			lonf = float(vals[2])
			tomof = float(vals[3])
			if (lat==latf):
				outfile.write('%g %g %g %g\n' %(radf,lonf,latf,tomof))

		print 'Sweep = %g' %i
		i+=1

	outfile.close()


#-----------------------------------------------------------
#Main code (needs to be organized better)
#-----------------------------------------------------------

#First, convert the tomography file to colatitude, and rewrite in the correct order
print '-----------------------------------------------------------'
print 'Converting to colatitude and reordering file for use with ASPECT'
print '-----------------------------------------------------------'
infile = open('DNA13_SVJtest.dat','r')
lines = infile.readlines()
infile.close()

#Determine the grid spacings 
tmpradius = []
for line in lines:
	vals = line.split('   ')
	R = float(vals[0])

	if R in tmpradius:
		break 
	else:
		tmpradius.append(R)

tmplat = []
for line in lines[0:len(lines):len(tmpradius)]:
	vals = line.split('   ')
	lat = 90-float(vals[1])

	if lat in tmplat:
		break 
	else:
		tmplat.append(lat)


tmplon = []
for line in lines[0:len(lines):len(tmpradius)*len(tmplat)]:
	vals = line.split('   ')
	lon = float(vals[2])

	if lon in tmplon:
		break 
	else:
		tmplon.append(lon)

colatconvert(tmplat,lines)

#-----------------------------------------------------------
#Read input tomography file and get information about the extent of the grid
#We convert immediately to radians, and colatitude (needed by ASPECT later)
#-----------------------------------------------------------
print '-----------------------------------------------------------'
print 'Generating interpolation grids'
print '-----------------------------------------------------------'
datafile=open('exout_colat_tmp.dat','r')
lines = datafile.readlines()
datafile.close()

#First, load the data
DATAVEC = np.zeros(len(lines))
i=0
for line in lines:
	vals = line.split(' ')
	VAL = float(vals[3])
	DATAVEC[i] = VAL
	i+= 1

#Determine the grid spacings 
tmpradius = []
for line in lines:
	vals = line.split(' ')
	R = float(vals[0])*1000.0

	if R in tmpradius:
		break 
	else:
		tmpradius.append(R)

tmplon = []
for line in lines[0:len(lines):len(tmpradius)]:
	vals = line.split(' ')
	lon = (float(vals[1]))*(np.pi/180.0)

	if lon in tmplon:
		break 
	else:
		tmplon.append(lon)


#This is colatitude
tmplat = []
for line in lines[0:len(lines):len(tmpradius)*len(tmplon)]:
	vals = line.split(' ')
	lat = float(vals[2])*(np.pi/180.0)

	if lat in tmplat:
		break 
	else:
		tmplat.append(lat)

#-----------------------------------------------------------
#Determine the number of points, minimum and maximium and the spacing
#-----------------------------------------------------------

#Determine the number of points in lat,lon and radius
Rpoints = len(tmpradius)
Latpoints = len(tmplat)
Lonpoints = len(tmplon)

Rpointinterval = abs(tmpradius[0]-tmpradius[1])
Latpointinterval = abs(tmplat[0]-tmplat[1])
Lonpointinterval = abs(tmplon[0]-tmplon[1])

Radiusbounds = [min(tmpradius),max(tmpradius)]
Latbounds = [min(tmplat),max(tmplat)]
Lonbounds = [min(tmplon),max(tmplon)]

#Generate the grids
RGRD = np.array(tmpradius)
LATGRD = np.array(tmplat)
LONGRD = np.array(tmplon)

#Data must be a 3D array
DATA = np.zeros((Rpoints,Lonpoints,Latpoints))

for i in range(Latpoints):
	for j in range(Lonpoints):
		for k in range(Rpoints):
			counter = i*Lonpoints*Rpoints + j*Rpoints + k
			DATA[k,j,i] = DATAVEC[counter]

print '---------------------------------------'
print 'Using the following grids'
print '---------------------------------------'

print 'LATITUDE'
print LATGRD
print 'RADIUS'
print RGRD
print 'LONGITUDE'
print LONGRD

interpfunc = RegularGridInterpolator((RGRD,LONGRD,LATGRD),DATA,method='nearest')

#-----------------------------------------------------------
#Create the global grid: Should be more sparsely spaced points than the regional tomo grid:
#whatever is needed to capture desired resolution in the global tomography
#-----------------------------------------------------------

minrad = 3490000.0
maxrad = 6370000.0

#parameters for the global grid
stepdivider = 2.0

radstepnumber = int(((maxrad-minrad)/Rpointinterval)/(0.5*stepdivider))
latstepnumber = int(((np.pi)/Latpointinterval)/(stepdivider))
lonstepnumber = int(((2*np.pi)/Lonpointinterval)/(stepdivider))

#global grid, in spherical coodinates

radgrid = np.linspace(minrad,maxrad,radstepnumber)

#this need to be colatitude!!!
latgrid = np.linspace(0,np.pi,latstepnumber)
longrid = np.linspace(0,2*np.pi,lonstepnumber)

outfileglobe = open('ASPECT_tomoregion_SVJTEMP_absolutetemps.dat','w')

tempdata = []
tempX = []
tempY = []
tempZ = []

#-----------------------------------------------------------
#Loop through the global tomography grid in 3D, and append points to cloud. Do not append points in the high resolution region: this comes later
#----------------------------------------------------------

#generate the meterial property interpolation functions, from Ian's curves 
f_temps,f_rhos,f_Ks,f_As,f_vs,fmeanvs,f_etas = AT.generateinterpfuncs()


#Remember, ASPECT requires points in radius, longitude, latitude; where all are ascending. Radius ascends first, then longitude, then latitude
outfileglobe.write('# POINTS: %i %i %i\n' %(radstepnumber,lonstepnumber,latstepnumber))
for lat in latgrid:
	for lon in longrid:
		for rad in radgrid:


  			### Here we've assumed that absolute velocity has been extracted: Find the material properties that correspond to this depth value
	  		tr = f_temps(rad)
	  		tA = f_As(rad)
	  		etaval = f_etas(rad) #get the dln(rho)/dln(V) value for this depth

	  		#Check if we're inside the domain 
	  		if (Radiusbounds[0] <= rad <= Radiusbounds[1]) and (Latbounds[0] <= lat <= Latbounds[1]) and (Lonbounds[0] <= lon <= Lonbounds[1]):

	  			if rad < 6270000: #Must be less than 100 km from the surface to extract from the interpolation function

					interpoint = interpfunc([rad,lon,lat]) #This should be some measure of velocity (a percentage perturbation, in this case)

					#meanvs = fmeanvs(rad) #The global model mean velocity at this depth
					
					#print 'Depth: %g' %(maxrad*1000-rad)
					#print 'model vel: %g' %tv
					
					#print 'interp vel: %g' %interpoint
					#print '-------------------'

					#Method 1: Simple: Don't use shear modulus to determine Rhodiff (seems to be more accurate!)

					#interpoint = 0 #Used if we just want the 3D radial temperature perturbation 
					temp = AT.ReturnTempPert(rad, interpoint, tr, tA, eta=etaval)

				else:
					Svel = 0.0
					temp = AT.ReturnTempPert(rad,Svel,tr,tA)


            #If we're not in the domain, just use the radial (adiabatic) temperature profile
			else:
				Svel = 0.0 #estimate of S velocity for the 'uppermost' part of the model - really the lithopshere since its anything above 100km. Set at 0.0 so that we get the radial
				#profile geotherm here
				temp = AT.ReturnTempPert(rad, Svel, tr, tA, eta=etaval)
				#interpoint = 5.0

  			outfileglobe.write('%g %g %g %g\n' %(rad,lon,lat,temp))

outfileglobe.close()

#Make a point cloud file containing the full grid
#X,Y,Z = spheretocart(np.array(tempX),np.array(tempY),np.array(tempZ))
#pointsToVTK("./ASPECT_TempSVJ_input",X,Y,Z,data = {"Pertb": np.array(tempdata)})


#-----------------------------------------------------------
#Now we loop though the points in the tomography datafile, appending to the point cloud to create a high-density region
#-----------------------------------------------------------

llen = len(lines)

#Make vectors for input into vtk program
Tomovals = np.zeros(llen)
Radii = np.zeros(llen)
Colat = np.zeros(llen)
Lon = np.zeros(llen)

i=0

for line in lines:
	vals = line.split(' ')
	R = float(vals[0])
	lat = float(vals[2])
	lon = float(vals[1])
	tomoval = float(vals[3]) #this is going to be in anomaly value

	#Write to ASPECT input
	r,rlon,rlat = converttoASPECT(R,lon,lat)
	Radii[i] = r
	Colat[i] = rlat
	Lon[i] = rlon
	Tomovals[i] = tomoval
	i+=1

Xnew,Ynew,Znew = spheretocart(Lon,Colat,Radii)

pointsToVTK("./TempTomoTestSVJ",Xnew,Ynew,Znew,data = {"Pertb": Tomovals})

