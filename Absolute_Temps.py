#!/usr/bin/env python

#This was decompiled from a .pyc file using
#uncompyler.py Absolute_Temps.pyc > recovered.py
#When the source code got writeen over

import numpy as np
from scipy import interpolate

def extract_vectors(inputfile):
    """Extract vectors from the model text file, ready to plot"""

    infile = open(inputfile, 'r')
    lines = infile.readlines()
    infile.close()
    rad = []
    grav = []
    pressure = []
    temp = []
    density = []
    vp = []
    vs = []
    K = []
    G = []
    alpha = []
    for line in lines[1:]:
        vals = line.split()
        rad.append(float(vals[0]))
        grav.append(float(vals[1]))
        pressure.append(float(vals[2]))
        temp.append(float(vals[3]))
        density.append(float(vals[4]))
        vp.append(float(vals[5]))
        vs.append(float(vals[6]))
        K.append(float(vals[7]))
        G.append(float(vals[8]))
        alpha.append(float(vals[9]))

    return (np.array(rad),
     np.array(temp),
     np.array(density),
     np.array(G),
     np.array(alpha),
     np.array(vs))



def extractmean(meansfile):
    """Extract the global mean of barbara's tomography - this is just a 1D profile, against which perurbations can be compared"""
    infile = open(meansfile, 'r')
    lines = infile.readlines()
    infile.close()
    Erad = 6371
    rads = []
    smeans = []
    for line in lines:
        vals = line.split()
        rad = (Erad - float(vals[0])) * 1000.0
        smean = float(vals[1])
        rads.append(rad)
        smeans.append(smean)

    return (np.array(list(reversed(rads))), np.array(list(reversed(smeans))))



def extracteta(etafile):
    """Extract variaiton in the parameter eta, which corresponds to dlin(rho)/dln(Vs) with depth. These values should come from 
    Karato et al, 1993"""
    infile = open(etafile, 'r')
    lines = infile.readlines()
    infile.close()
    Erad = 6371
    rads = []
    etas = []
    for line in lines:
        vals = line.split()
        rad = (Erad - float(vals[0])) * 1000.0
        eta = float(vals[1])
        rads.append(rad)
        etas.append(eta)

    print rads,
    print etas
    return (np.array(list(reversed(rads))), np.array(list(reversed(etas))))



def generateinterpfuncs():
    """
            Generate 1D interpolation functions from the input data
    
            We need to determine density from velocity, which requires knowlege of K, the the shear modulus
            We then map from density to temperature pertubation, which needs knowlege of alpha, T0, rho and Rho0 - all of which we get from the profiles
            """
    (R, T, RHO, G, A, VS,) = extract_vectors('/Users/rmartinshort/Documents/Berkeley/ASPECT/Tomo_data_extract/Global_tomo/Extract_from_SEMUCBWM1/Material_profiles/mantle_model.txt')
    (Rmeans, Smeans,) = extractmean('/Users/rmartinshort/Documents/Berkeley/ASPECT/Tomo_data_extract/Global_tomo/Extract_from_SEMUCBWM1/Barb_Tomo_means.dat')
    (Retas, etas,) = extracteta('/Users/rmartinshort/Documents/Berkeley/ASPECT/Tomo_data_extract/Global_tomo/Extract_from_SEMUCBWM1/Material_profiles/Rho_vs_scaling.dat')
    f_temps = interpolate.interp1d(R, T)
    f_rhos = interpolate.interp1d(R, RHO)
    f_Ks = interpolate.interp1d(R, G)
    f_As = interpolate.interp1d(R, A)
    f_vs = interpolate.interp1d(R, VS)
    fmeanvs = interpolate.interp1d(Rmeans, Smeans)
    f_etas = interpolate.interp1d(Retas, etas)

    return (f_temps,
     f_rhos,
     f_Ks,
     f_As,
     f_vs,
     fmeanvs,
     f_etas)



def ReturnTempPert(Rad, Spertin, Tr, Ar, eta = 0.25):
    """Takes a pertubation in S velocity value and used interpolation functions to determine the temperature at that depth
    The value of eta = 0.25 is suggested by the ASPECT manual, but the interpolated value from Karato (1993) could be used instead"""

    drho = eta * (float(Spertin) / 100.0)
    Tpert = Tr - 1 / Ar * drho

    if Tpert > 2900:
        Tpert = 2900

    return Tpert



def ReturnTempAbsolute(Rad, Svelin, Smeanvelin, Tr, Kr, Ar):
    """Takes an absolute S velocity value and uses the interpolation functions to determine the temperature at that depth. This method uses information
    about the shear modulus of the material with depth to determin the density difference, so it should be more accurate than the above"""

    rhopert = Kr / (Svelin * 1000) ** 2
    rhopertmean = Kr / (Smeanvelin * 1000) ** 2
    drho = rhopertmean - rhopert
    Tpert = Tr + drho / (rhopert * Ar)


    if Tpert > 2900:
        Tpert = 2900

    return Tpert



# decompiled 1 files: 1 okay, 0 failed, 0 verify failed
# 2016.02.26 16:15:05 PST
