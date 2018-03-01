import os, matgen as mg

A = mg.getSemiconductorsWithBS()
mpids = [Ai['material_id'] for Ai in A]

# stability
eh = {Ai['material_id']: Ai['e_above_hull'] for Ai in A}

# carrier masses
me_max = {}
me_min = {}
me_mean = {}
mh_max = {}
mh_min = {}
mh_mean = {}
for i, mpid in enumerate(mpids):
    print('***** {} / {} *****'.format(i+1, len(mpids)))
    try:
        d = mg.jsonLoad(os.path.join(mg.ad, 'MASS', mpid))
        me_min[mpid] = d['me_min']
        me_max[mpid] = d['me_max']
        me_mean[mpid] = d['me_mean']
        mh_min[mpid] = d['mh_min']
        mh_max[mpid] = d['mh_max']
        mh_mean[mpid] = d['mh_mean']
    except:
        print('FAILED FOR {}'.format(mpid))

# band gap info
gt = {}
eg = {}
dg = {}
for i, mpid in enumerate(mpids):
	print('***** {} / {} *****'.format(i+1, len(mpids)))
	try:
		d = mg.getBS(mpid)
		gt[mpid] = d.get_band_gap()['direct']
		eg[mpid] = d.get_band_gap()['energy']
		dg[mpid] = d.get_direct_band_gap()
	except:
		print('FAILED FOR {}'.format(mpid))

# combine into a single object for easy storage
data = {'eh': eh, 'gt': gt, 'dg': dg, 'eg': eg, 'me_max': me_max,
	'me_min': me_min, 'me_mean': me_mean, 'mh_max': mh_max, 'mh_min': mh_min, 'mh_mean': mh_mean}

# find common entries (keys)
K = set.intersection(*[set(list(p.keys())) for p in data.values()])

