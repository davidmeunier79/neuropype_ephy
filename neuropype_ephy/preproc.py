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

def preprocess_ts(ts_file,orig_channel_names_file,orig_channel_coords_file,orig_sfreq, down_sfreq,prefiltered = False):
    
	from mne.io import RawArray	
	
	from mne import create_info
	
	import os
	import numpy as np


	


        
	#### load electrode names
	elec_names = [line.strip() for line in open(orig_channel_names_file)]
	#print elec_names

	### save electrode locations
	elec_loc = np.loadtxt(orig_channel_coords_file)
	#print elec_loc

        ### no modification on electrode names and locations
        correct_elec_loc = elec_loc
        correct_elec_names = elec_names

        print len(correct_elec_names)
        print len(correct_elec_loc)
        
        
        ### save electrode locations	
	channel_coords_file = os.path.abspath("correct_channel_coords.txt")
	np.savetxt(channel_coords_file ,correct_elec_loc , fmt = '%s')

	#### save electrode names
	channel_names_file = os.path.abspath("correct_channel_names.txt")
	np.savetxt(channel_names_file,correct_elec_names , fmt = '%s')

        

        ##### downsampling on data
        ts = np.load(ts_file)
        
        print ts.shape
        
        0/0
        
        
        raw = RawArray(ts, info = create_info(ch_names = elec_names, sfreq = orig_sfreq))
        
        indexes_good_elec = np.arange(len(elec_names))
        
        print indexes_good_elec
        
        if prefiltered == False:
		raw.filter(l_freq = None, h_freq = down_sfreq, picks = indexes_good_elec)

	raw.resample(sfreq = down_sfreq,npad = 100)
	
	downsampled_ts,times = raw[:,:]


        print downsampled_ts.shape
        
        
	downsampled_ts_file = os.path.abspath("downsampled_ts.npy")

	np.save(downsampled_ts_file,downsampled_ts)

        print raw.info['sfreq']
        
	return downsampled_ts_file,channel_coords_file,channel_names_file,raw.info['sfreq']
