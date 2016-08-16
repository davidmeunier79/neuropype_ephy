def compute_and_save_psd(epochs_fname,  fmin=0, fmax=120, method='welch',
						 n_fft=256, n_overlap=0, picks=None,
						 proj=False, n_jobs=1, verbose=None):
	"""
	Load epochs from file, 
	compute psd and save the result in numpy arrays
	"""
	import numpy as np
	import os
	from mne import read_epochs
	epochs = read_epochs(epochs_fname)
	if method == 'welch':
		from mne.time_frequency import psd_welch
		psds, freqs = psd_welch(epochs)
	elif method == 'multitaper':
		from mne.time_frequency import psd_multitaper
		psds, freqs = psd_multitaper(epochs)
	else:
		raise Exception('nonexistent method for psd computation')
	path, name = os.path.split(epochs_fname)
	base, ext = os.path.splitext(name)
	psds_fname = base + '-psds.npy'
	freqs_fname = base + '-freqs.npy'
	psds_file = os.path.abspath(psds_fname)
	freqs_file = os.path.abspath(freqs_fname)
	np.save(psds_file, psds)
	np.save(freqs_file, freqs)
	return psds_file, freqs_file
