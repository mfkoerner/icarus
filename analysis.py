

import os
from electrons import np, plt
import photons as ph
from matgen import td, ad

bd = os.path.join(td, 'General-Direct')

def Jmax(mpid):
	''' '''
	hnu, alpha, eps, N, Eg = ph.optical(os.path.join(bd, mpid, 'absorb'))
	xx, yySun, yyMat, ss = ph.FOM(hnu, alpha, Eg)
	return ss[-1]/1e8, Eg