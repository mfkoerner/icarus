# IMPORT PACKAGES
import os, importlib as imp, matplotlib.gridspec as gs, pandas as pd
try:
	import photons as ph
except:					#This fails whenever the top one fails
	imp.reload(ph)
plt = ph.plt
np = ph.np
import matgen as mg

# USER SETTINGS
Egmin = 0.
xl = (200, 3000)
"""
topList explanation:
this should have a specific file structure of compoundname/absorb/{DOSCAR,IBZKPT,OUTCAR} at the minimum
don't list the compoundnames, just the dirctory holdin those. Most likely you only want one entry in topList
"""
topList = [
"file_directories_here"
	]
bdList, lblList = [], []
for batch in topList:
	dirList = next(os.walk(batch))[1]
	for compound in dirList:
		if 'absorb' in os.listdir(os.path.join(batch, compound)):
			bdList.append(os.path.join(batch, compound, 'absorb'))
			lblList.append(compound)

# DO STUFF
lamSun, ySun = ph.getSolarSpectrum()				#Solar spectrum numbers (irradiance vs wavelength)
fig = plt.figure()
grid = gs.GridSpec(3, 1, hspace=0.02, height_ratios=(1, 3, 3))	#Create grid for axes
ax0 = plt.subplot(grid[0])							#Axis0 for solar spectrum
ax0.plot(lamSun, ySun)
ax0.set_ylabel('irradiance\n(W m$^{-2}$ nm$^{-1}$)')
ax1 = plt.subplot(grid[1])							#Axis1 for absorption vs wavelength
ax2 = plt.subplot(grid[2])							#Axis2 for FOM
allEg = {}											#allEG is for storing Energy vector for each materials along dielectric function
allFOM = {}											#allFOM stores the FOM data
for bd, lbl in zip(bdList, lblList):				#bd is absorption directory (basedir), lbl is compound name
	# try:
	hnu, alpha, eps, N, Eg = ph.optical(bd)			#energy, absorption coefficient, effective eigenvalues of dielectric matrix, dielectric constant, band gap
	if Eg>Egmin: 									#Currently just requires positive band gap
		ind = hnu>Eg 								#Matrix of TRUE for frequencies with energy above the band gap, and FALSE for energies below
		xx, yySun, yyMat, ss = ph.FOM(hnu, alpha, Eg)	#Wavelength, power/(area*wavelength), absorption (by wavelength), total absorbed light (cumulatively integrated)
		lh = ph.plotAbsorption(ax1, hnu[ind], alpha[ind], wavelength=True, xl=(100,4000))	#
		if 'mp-' in lbl:
			lh.set_label('{} ({})'.format(mg.getFormula(lbl), lbl))
		else:
			lh.set_label(lbl)
		print(lbl, Eg, ss[-1]/1e8, len(hnu))
		allEg[lbl] = Eg
		allFOM[lbl] = ss[-1]
		if len(hnu) < 500:
			lh.set_color((0.5,)*3)
			lh.set_linestyle('--')
			lh2 = ax2.plot(xx[1:], ss/1e8, '--', c=(0.5,)*3, label=lbl)
		else:
			lh2 = ax2.plot(xx[1:], ss/1e8, label=lbl)
	else:
		print('band gap for {} too small (under {})'.format(lbl,Egmin))
	# except:
		# print('FAILED TO READ {}'.format(bd))
ax0.set_xticklabels('')
ax1.set_xticklabels('')
ax1.set_xlabel('')
ax2.set_xlabel('$\\lambda$ (nm)')
ax2.set_ylabel('FOM / 10$^8$')
ax0.set_xlim(xl)
ax0.set_ylim(0, 2)
ax1.set_xlim(xl)
ax1.set_ylim(1e3, 1e7)
ax2.set_xlim(xl)
ax2.set_ylim(0, 6)
fig.show()



fig2, ax = plt.subplots()
K = allEg.keys()
Eg = np.array([allEg[ki] for ki in K])
FOM = np.array([allFOM[ki]/1e8 for ki in K])
combined = {ki: allEg[ki]*allFOM[ki]/1e8 for ki in K}

d1 = pd.read_csv('data/autoparse/Eff-Mass/MASS.txt', index_col=0)
d2 = pd.read_csv('data/autoparse/Gap-and-Type/GAP.txt', index_col=0)
d2 = d2.drop('formula', axis=1)
d3 = pd.DataFrame.from_dict(combined, orient='index')
d3 = d3.rename(columns={0: 'Eg*FOM'})
d = d1.merge(d2, how='inner', left_index=True, right_index=True)
d = d.merge(d3, how='inner', left_index=True, right_index=True)

# ax.plot(Eg, FOM, 'o', label='FOM / 10$^8$')
# ax.plot(Eg, [combined[ki] for ki in K], '^', label='$E_g$ * FOM / 10$^8$')

ax.plot(d['me_mean'], d['Eg*FOM'], 'o', label='$m_e$')
ax.plot(-d['mh_mean'], d['Eg*FOM'], '^', label='$m_h$')
ax.set_xlim(0, 1)
ax.legend()
fig2.show()
input('Press <ENTER> to exit')
