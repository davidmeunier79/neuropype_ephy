# -*- coding: utf-8 -*-

def preprocess_fif_to_ts(fif_file):

	import os
	import numpy as np

	from mne.io import Raw	
	#from mne.io import RawFIF	## was working in previous versions ofpyMNE
	from nipype.utils.filemanip import split_filename as split_f

	subj_path,basename,ext = split_f(fif_file)

	print fif_file
	
	raw = Raw(fif_file,preload = True)

	print raw

	print len(raw.ch_names)

	select_sensors, = np.where(np.array([ch_name[0] == 'M' for ch_name in raw.ch_names],dtype = 'bool') == True)

	### save electrode locations
	sens_loc = [raw.info['chs'][i]['loc'][:3] for i in select_sensors]
	sens_loc = np.array(sens_loc)

	channel_coords_file = os.path.abspath("correct_channel_coords.txt")
	np.savetxt(channel_coords_file ,sens_loc , fmt = '%s')

	print sens_loc

	### save electrode names
	sens_names = np.array([raw.ch_names[pos] for pos in select_sensors],dtype = "str")

	channel_names_file = os.path.abspath("correct_channel_names.txt")
	np.savetxt(channel_names_file,sens_names , fmt = '%s')

	### filtering + downsampling

	data,times = raw[select_sensors,:]
	print data.shape
	print raw.info['sfreq']

	raw.filter(l_freq = None, h_freq = 300,picks = select_sensors)

	raw.resample(sfreq = 300,npad = 0,stim_picks = select_sensors)


	### save data
	data,times = raw[select_sensors,:]
	print data.shape
	print raw.info['sfreq']
	#0/0


	ts_file = os.path.abspath(basename +".npy")

	np.save(ts_file,data)

	return ts_file,channel_coords_file,channel_names_file,raw.info['sfreq']
