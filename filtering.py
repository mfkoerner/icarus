import matgen as mg

def filter(data, remMetals=True, spinPol=False, directOnly=True, ehmax=0.05, dgDeltaMax=None, memax=None, mhmax=None):
	''' '''
	K = list(set.intersection(*[set(list(p.keys())) for p in data.values()]))
	print('{} COMMON ENTRIES FOUND ACROSS ALL DATA TYPES'.format(len(K)))
	# remove metals (some erroneously included from database query for semiconductors)
	if remMetals:
		K = [k for k in K if data['eg'][k] != 0]
	# eliminate indirect gap materials?
	if directOnly:
		K = [k for k in K if data['gt'][k]]
	# crude stability criterion
	K = [k for k in K if data['eh'][k] < ehmax]
	# direct transition delta criterion
	if dgDeltaMax is not None:
		K = [k for k in K if data['dg'][k] < data['eg'][k] + dgDeltaMax]
	# effective mass criteria
	if memax:
		K = [k for k in K if data['me_mean'][k] < memax]
	if mhmax:
		K = [k for k in K if data['mh_mean'][k] < mhmax]
	# remove magnetic materials
	if not spinPol:
		K = [k for k in K if not mg.getBS(k).is_spin_polarized]
	print('{} ENTRIES REMAIN AFTER PROPERTY FILTERING'.format(len(K)))
	return set(K)
