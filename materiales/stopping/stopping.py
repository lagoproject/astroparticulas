# -*- encoding: utf-8 -*-
import sys
import math

'''
Very simple code to perform numerical integrations of stopping power curves (obtained from, eg, NIST) to calculate total ranges as a function of particle initial momentum for different materials or densities, under the CSD approximation. It produces also bragg curves for each injected particle. 

Original version and further modifications written in perl by Hernán Asorey, Centro Atómico Bariloche, (C) 2012, asoreyh@gmail.com.

New version adapted and improvements by Marta Valencia, Universidad Industrial de Santander, (C) 2015.
'''

# parameters 
colp    = 0 	# table column for kinetical energy (T in MeV)
cold    = 7 	# table column for de/dx in MeV cm2 / g
m       = 105.7	# particle mass, in MeV, needed for calculate the total energy, 
                # E = sqrt(p2+m2), and T=E-m
l       = 0.1 	# step in cm, same units given in tables by NIST. Usually, 1 mm for small objects
density = 0. 	# density parameter, g/cm3, if 0 it will try to find density at NIST file
pi      = 5000. # current particle momentum (initial), in MeV
pimax   = 10000. # max particle momentum, should be > pi
pistep  = 0.05  # step for particle momentum, in fraction, ie, 0.05 means momentum 
                # will be increased in steps of 5%, pi *= 1+pistep, 
                # sarting from pi and ending at pimax
# files
comment = "#"	                        # comment identifier
filein  = "muon_water_dedx.dat"         # dE/dX ASCII file
fileout = "muon_water"                  # base filename for results
bragg   = True                          # if true produce bragg curves for each initial momentum. warning, could produce several files

# arrays
dk = [] # kinetical energy from table, in MeV
de = [] # stopping power from table, in MeV cm2 / g
d  = [] # auxiliar array

# let's work
data = open(filein, "r")
filerange = '{0}_ranges.dat' .format(fileout)
f = open(filerange, "w")


# reading ASCII table
for readline in data:
    d = readline.split()
    if not density: 
        if "density" in readline:
            density = d[d.index("density") + 2]
    if d[0] != comment:
        if density == 0.:
            f.write("\nERROR: Density parameter is needed.\n")
            sys.exit()  # do not find density or/and density was not defined...
        # everything fine, just fill the arrays
        dk.append(d[colp])
        c = float(d[cold])*float(density) # assuming constant density
        de.append(float(c))

# starting output file
f.write("# Momemtum(MeV/c) Total_Energy(MeV) Kinetic_Energy(MeV) Range(cm)\n")

while (pi < pimax):
    rge = 0.
    t = len(dk) - 1 # energy index
    e = math.sqrt(pi**2 + m**2) #from p to energy
    k = e - m
    f.write("%f %f %f " % (pi, e, k))
    if bragg:
        filebragg = '{0}_bragg_{1}.dat' .format(fileout, pi)
        g = open(filebragg, "w")
        g.write("# Current_Energy(MeV) Current_K_energy(MeV) dedx(MeV_cm2/g) Distance(cm)\n")
    while (e > m): # moving the particle if k>0
        # 1. find $k in dedx data
        while (t >= 0 and float(dk[t]) > k):
            t -= 1
        if (t < 0):
            t = 0
        # 2. interpolate to get dedx
        x1 = float(dk[t])
        x2 = float(dk[t+1])
        y1 = float(de[t])
        y2 = float(de[t+1])
        m1 = (y2 - y1) / (x2 - x1)
        m0 = y1 - m1 * x1
        # 3. for this energy, de/dx is
        dedx = (m1 * k + m0) * l
        if dedx < 0.: # that's it, typical low energy problems...
            e=m
            break
        # 4. and the range up to here is
        rge += l
        # 5. print results
        if bragg:
            g.write("%f %f %f %f\n" % (e, k, dedx, rge))
        # 6. then, the new energies are
        e -= dedx
        if (e < m):  # this is the end, my friend
            e = m
        k = e - m
        # end of particle loop
    
    f.write ("%f\n" % (rge))
    if bragg:
        g.close()
    pi *= (1. + pistep)
# end of momentum loop
f.close()
