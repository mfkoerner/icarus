########################################
# written for Python 3                 #
# by Doug Fabini (fabini@mrl.ucsb.edu) #
########################################

'''

 This script requires the following files to be located in 'baseDir':
 - IBZKPT (to extract number of k points) POSSIBLY NO LONGER NEEDED
 - DOSCAR (to extract bandgap)
 - OUTCAR (to extract dielectric properties and energy resolution)

Currently only handles an isotropic equivalent for the dielectric / absorption tensors.

'''



# import packages, apply stylesheet
import config
import os
import numpy as np
import matplotlib.pyplot as plt
import electrons as el

plt.style.use(config.plot_style)


# ****************** #
# DATA I/O FUNCTIONS #
# ****************** #

def getNkPts(bd):
	''' Parse OUTCAR for number of k-points '''
	fname = os.path.join(bd, 'OUTCAR')
	# print(fname) #debug line
	with open(fname, 'r') as f:
		for line in f:
			if 'irreducible k-points' in line:
				# print(line)    #debug line
				return int(line.split()[1])
				break

def getDielectric(bd):
	''' Parse OUTCAR for dielectric properties, convert to appropriate form '''
	fname = os.path.join(bd, 'OUTCAR')
	with open(fname, 'r') as f:
		raw = []
		lnImag, lnReal = 0, 0
		for i, line in enumerate(f):
			raw.append(line)
			if 'NEDOS' in line: 					#This an below find number points per section and start of lines for sections
				NEDOS = int(line.split()[5])
			if 'IMAGINARY DIELECTRIC' in line and lnImag is 0:	#Selecting the first set of Dielectric numbers from VASP
				lnImag = i
			if 'REAL DIELECTRIC' in line and lnReal is 0:
				lnReal = i
	EepsRe, EepsIm = [], []
	for i in range(lnImag+3,lnImag+NEDOS+3):		#All of the imaginary dielectric components (NEDOS components and start point of lnImag+3)
		if len(raw[i]) < 5:							#Checking for early termination of DIELECTRIC DATA (printing to output)
			print('DIELECTRIC DATA TERMINATED AT ONLY {} POINTS'.format(i-lnImag-3))
			break
		EepsIm.append([float(ri) for ri in raw[i].strip('\n').split()])	#Energy (frequency) then X,Y,Z,XY,YZ,ZX for imaginary component
	E = np.array(EepsIm)[:,0]						#Energies pulled from first part of EepsIm
	epsIm = np.array([isotropic(row[1:]) for row in EepsIm])	#epsIm is the isotropic equivilent values for each energy
	for i in range(lnReal+3,lnReal+NEDOS+3):
		if len(raw[i]) < 5:
			# print('DIELECTRIC DATA TERMINATED AT ONLY {} POINTS'.format(i-lnReal-3))
			break
		EepsRe.append([float(ri) for ri in raw[i].strip('\n').split()])	#Real part from above
	epsRe = np.array([isotropic(row[1:]) for row in EepsRe])	#Real part for epsIm, this time is epsRe
	return E, epsRe + 1j*epsIm 					#Returns list of isotropic equivalent values

def saveResults(bd, E, alpha, eps):
	''' Store absorption coefficient and dielectric function '''
	out = np.hstack((E, alpha, eps.real, eps.imag))
	out = np.reshape(out, (-1, 4), order='F')
	np.savetxt(os.path.join(bd, 'optical.csv'), out, header='h*nu (eV), alpha_iso (cm^-1), Re[eps_iso] (eps_0), Im[eps_iso] (eps_0)')

def getSolarSpectrum():
	''' Get direct+diffuse solar irradiance at global tilt, ASTM G173-03 '''
	d = np.loadtxt('data/ASTMG173.dat')
	return d[:,0], d[:,2]


# ****************** #
# ANALYSIS FUNCTIONS #
# ****************** #

def nm2eV(lam):
	''' Convert wavelength in nm to energy in eV '''
	h = 4.136e-15 # Planck constant, eV / s
	c = 2.998e8   # speed of light, m / s
	return h*c/(lam*1e-9)

def eV2nm(hnu):
	''' Convert energy in eV to wavelength in nm '''
	h = 4.136e-15 # Planck constant, eV / s
	c = 2.998e8   # speed of light, m / s
	return h*c/hnu*1e9

def isotropic(sixElements):
	''' Returns an isotropic equivalent value for a symmetric 3x3 matrix '''
	xx, yy, zz, xy, yz, zx = sixElements
	A = np.array([[xx, xy, zx], [xy, yy, yz], [zx, yz, zz]])
	eigval, _ = np.linalg.eigh(A)
	return np.mean(eigval)

def dielec2optical(hnu, eps):
	''' Calculate complex refractive index and absorption coefficient from dielectric function '''
	h = 4.136e-15 # Planck constant, eV / s
	c = 2.998e8   # speed of light, m / s
	N = np.sqrt(eps)
	alpha = 4*np.pi/(h*c)*hnu*N.imag/100 # divisor of 100 takes from m-1 to cm-1
	return N, alpha

def FOM(hnu, alpha, Eg):

	xx = np.linspace(100, eV2nm(Eg), int(1e4))  							#proper range of light to think about (100 nm [13eV] to band gap wavelength)
	xSun, ySun = getSolarSpectrum() 										#xSun -> wavelength of sun, ySun -> intensity of sun
	yySun = np.interp(xx, xSun, ySun) 										#ySun calculated at the points for xx (so that we have the right resolution)
	yyMat = np.interp(xx, np.flipud(eV2nm(hnu[1:])), np.flipud(alpha[1:])) 	#absorption as a function of wavelength
	from scipy.integrate import cumtrapz 									#Trapezoidal numeric integration
	return xx, yySun, yyMat, cumtrapz(yySun*yyMat, xx) 						#FOM is the last value, which is integral of sum intensity time absorption along wavel



# ****************** #
# PLOTTING FUNCTIONS #
# ****************** #

def plotDielectric(ax, E, eps, N, El=(0, 10)):
	''' Plot complex dielectric function and complex refractive index '''
	ax.plot(E, eps.real, label='$\\epsilon_r\\prime$')
	ax.plot(E, eps.imag, label='$\\epsilon_r\\prime\\prime$')
	ax.plot(E, N.real, label='$n$')
	ax.plot(E, N.imag, label='$k$')
	ax.set_xlim(El)
	ax.set_xlabel('$h\\nu$ (eV)')
	ax.legend()

def plotAbsorption(ax, hnu, alpha, xl=(0, 4), yl=(1e2, 1e7), rel2eg=None, lbl=None, wavelength=False):
	''' Plot absorption coefficient '''
	if wavelength:
		if rel2eg is not None:
			raise Exception('Relative to gap option not available when plotting by wavelength')
		lh, = ax.semilogy(eV2nm(hnu), alpha, '.-', label=lbl)
		ax.set_xlabel('$\\lambda$ (nm)')
	elif not wavelength and rel2eg is None:
		lh, = ax.semilogy(hnu, alpha, '.-', label=lbl)
		ax.set_xlabel('$h\\nu$ (eV)')
	else:
		lh, = ax.semilogy(hnu-rel2eg, alpha, '.-', label=lbl)
		ax.set_xlabel('$h\\nu-E_g$ (eV)')
	ax.set_xlim(xl)
	ax.set_ylim(yl)
	ax.set_ylabel('$\\alpha$ (cm$^{-1}$)')
	return lh



# ********** #
# HIGH LEVEL #
# ********** #

def optical(bd):
	''' DESCRIPTION GOES HERE '''
	Nk = getNkPts(bd) 							#Gets number of irreducible kpoints but never uses it :O
	E, eps = getDielectric(bd)   				#Gets lists of E and equivilent eigenvalues (real + i*imag) for dialectric function
	N, alpha = dielec2optical(E, eps)			#N (dielectric constant) and alpha (absorption coefficient) from dielectric equivilent eigenvalues

	Edos, tdos = el.getTotalDOS(bd)				#arrays of len NEDOS with energy and DOS at that energy
	Eg = el.bandgap(Edos, tdos)					#Calculates bandgap from DOS data

	saveResults(bd, E, alpha, eps)				#Saves Energy, absorption, eigenvalue to basedir/optical.csv
	return E, alpha, eps, N, Eg 				#Returns Energy, absorption, eigenvalue, refractive index, bandgap
