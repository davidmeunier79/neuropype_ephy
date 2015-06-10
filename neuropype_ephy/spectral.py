# -*- coding: utf-8 -*-

import numpy as np

def spectral_proc(ts_file,sfreq,freq_band,freq_band_name):

	import numpy as np
	import os

	from params import con_method
	from mne.connectivity import spectral_connectivity

	data = np.load(ts_file)

	print data.shape
	print sfreq
	print freq_band
	print freq_band_name

	if len(data.shape) == 2:
		data = data.reshape(1,data.shape[0],data.shape[1])

	con_matrix, freqs, times, n_epochs, n_tapers  = spectral_connectivity(data, method=con_method, mode='multitaper', sfreq=sfreq, fmin= freq_band[0], fmax=freq_band[1], faverage=True, tmin=None,    mt_adaptive=False, n_jobs=1)

	con_matrix = np.array(con_matrix[:,:,0])

	print con_matrix.shape
	print np.min(con_matrix),np.max(con_matrix)
	
	conmat_file = os.path.abspath("conmat_"+ con_method + ".npy")

	np.save(conmat_file,con_matrix)

	return conmat_file


def multiple_spectral_proc(ts_file,sfreq,freq_band_name,freq_band,con_method):

	import numpy as np
	import os

	from params import con_method
	from mne.connectivity import spectral_connectivity

	all_data = np.load(ts_file)

	print all_data.shape
	
	#print sfreq
              
	print freq_band
	print freq_band_name
        
        if len(all_data.shape) != 3:
            print "Warning, all_data should have several samples"
            
            return []
        
        conmat_files = []
        
	for i in range(all_data.shape[0]):
	
                cur_data = all_data[i,:,:]
	
		data = cur_data.reshape(1,cur_data.shape[0],cur_data.shape[1])

                print data.shape
                
                con_matrix, freqs, times, n_epochs, n_tapers  = spectral_connectivity(data, method=con_method, mode='multitaper', sfreq=sfreq, fmin= freq_band[0], fmax=freq_band[1], faverage=True, tmin=None,    mt_adaptive=False, n_jobs=1)

                con_matrix = np.array(con_matrix[:,:,0])

                print con_matrix.shape
                print np.min(con_matrix),np.max(con_matrix)
                
                conmat_file = os.path.abspath("conmat_"+ con_method + "_" + str(i) + ".npy")

                np.save(conmat_file,con_matrix)

                conmat_files.append(conmat_file)
                
	return conmat_files

def epoched_multiple_spectral_proc(ts_file,sfreq,freq_band_name,freq_band,con_method,epoch_window_length):

	import numpy as np
	import os

	from mne.connectivity import spectral_connectivity

	all_data = np.load(ts_file)

	print all_data.shape
	
	#print sfreq
              
	print freq_band
	print freq_band_name
        
        if len(all_data.shape) != 3:
            print "Warning, all_data should have several samples"
            
            return []
        
        conmat_files = []
        
	for i in range(all_data.shape[0]):
	
                cur_data = all_data[i,:,:]
	
                print cur_data.shape
                    



                if epoch_window_length == None :
                    
                    data = cur_data.reshape(1,cur_data.shape[0],cur_data.shape[1])

                else: 
                        
                    nb_splits = cur_data.shape[1] // (epoch_window_length * sfreq)
                    
                    print "epoching data with {}s by window, resulting in {} epochs".format(epoch_window_length,nb_splits)
                    
                    list_epoched_data = np.array_split(cur_data,nb_splits,axis = 1)
                    
                    print len(list_epoched_data)
                    
                    data = np.array(list_epoched_data)
                    
                    print data.shape

                con_matrix, freqs, times, n_epochs, n_tapers  = spectral_connectivity(data, method=con_method, mode='multitaper', sfreq=sfreq, fmin= freq_band[0], fmax=freq_band[1], faverage=True, tmin=None,    mt_adaptive=False, n_jobs=1)

                print con_matrix.shape
                con_matrix = np.array(con_matrix[:,:,0])

                print con_matrix.shape
                print np.min(con_matrix),np.max(con_matrix)
                
                conmat_file = os.path.abspath("conmat_"+ con_method + "_" + str(i) + ".npy")

                np.save(conmat_file,con_matrix)

                conmat_files.append(conmat_file)
                
	return conmat_files





def compute_and_save_coherency_spectral_connectivity(data,con_method,sfreq,fmin,fmax,index = 0):
        
    import sys,os
    from mne.connectivity import spectral_connectivity
    
    import numpy as np

    if con_method in ['coh','cohy','imcoh']:
            
        print data.shape
        
        0/0
        
        if len(data.shape) < 3:
        	
        	data = data.reshape(1,data.shape[0],data.shape[1])
            
        con_matrix, freqs, times, n_epochs, n_tapers  = spectral_connectivity(data, method=con_method, mode='multitaper', sfreq=sfreq, fmin= fmin, fmax=fmax, faverage=True, tmin=None,    mt_adaptive=False, n_jobs=1)

        con_matrix = np.array(con_matrix[:,:,0])

        print con_matrix.shape
        print np.min(con_matrix),np.max(con_matrix)
        
        conmat_file = os.path.abspath("conmat" + str(index) + "_" + con_method + ".npy")
	
	np.save(conmat_file,con_matrix)

        return conmat_file
    
    else:
        
        print "warning, only work with coherency-based metrics"
        sys.exit()
        
def compute_and_save_phase_spectral_connectivity(epoched_data,con_method,sfreq,fmin,fmax,index = 0):
    
    import sys,os
    from mne.connectivity import spectral_connectivity
    
    import numpy as np

    if con_method in ['pli','plv','ppc' ,'pli','pli2_unbiased' ,'wpli' ,'wpli2_debiased']:
            
        if len(epoched_data.shape) < 3:
            print "warning, only work with epoched time series"
            sys.exit()
            
        con_matrix, freqs, times, n_epochs, n_tapers  = spectral_connectivity(epoched_data, method=con_method, mode='multitaper', sfreq=sfreq, fmin= fmin, fmax=fmax, faverage=True, tmin=None,    mt_adaptive=False, n_jobs=1)

        print con_matrix.shape
        
        con_matrix = np.array(con_matrix[:,:,0])

        print con_matrix.shape
        print np.min(con_matrix),np.max(con_matrix)
        
        conmat_file = os.path.abspath("conmat" + str(index) + "_" + con_method + ".npy")

        np.save(conmat_file,con_matrix)

        return conmat_file
    
    else:
        
        print "warning, only work with coherency-based metrics"
        sys.exit()
    
def epoched_spectral_proc(ts_file,sfreq,freq_band,freq_band_name,con_method,epoch_window_length):

	import numpy as np

	from neuropype_ephy.spectral import compute_and_save_coherency_spectral_connectivity,compute_and_save_phase_spectral_connectivity

	data = np.load(ts_file)

	print data.shape
	print sfreq
	print freq_band
	print freq_band_name

        if con_method in ['pli','plv','ppc' ,'pli','pli2_unbiased' ,'wpli' ,'wpli2_debiased']:
            
            if epoch_window_length == None:
                
                print "WARNING, phase-based metric will not work if epoch_window_length is not defined" 
            
            else:
                    
                nb_splits = data.shape[1] // (epoch_window_length * sfreq)
                
                print "epoching data with {}s by window, resulting in {} epochs".format(epoch_window_length,nb_splits)
                
                epoched_data = np.array(np.array_split(data,nb_splits,axis = 1))
                
                print epoched_data.shape

                conmat_file = compute_and_save_phase_spectral_connectivity(epoched_data=epoched_data, con_method=con_method, sfreq=sfreq, fmin= freq_band[0], fmax=freq_band[1])

                return conmat_file

        elif con_method in ['coh','cohy','imcoh']:
                
            if epoch_window_length == None:
                
                
                conmat_file = compute_and_save_coherency_spectral_connectivity(data=data,con_method=con_method,sfreq=sfreq,fmin = freq_band[0],fmax = freq_band[1])
                    
            
            else:
                
                print "Shape before splits:"
                print data.shape
                
                print  "sfreq:"
                print sfreq
                
                nb_splits = data.shape[1] // (epoch_window_length * sfreq)
                
                print "nb_splits:"
                print nb_splits
                
                print "epoching data with {}s by window, resulting in {} epochs".format(epoch_window_length,nb_splits)
                
                0/0
                
                
                epoched_data = np.array(np.array_split(data,nb_splits,axis = 1))
                
                print epoched_data.shape

                conmat_file = compute_and_save_coherency_spectral_connectivity(data=epoched_data, con_method=con_method, sfreq=sfreq, fmin= freq_band[0], fmax=freq_band[1])

                return conmat_file
                
                
                ### previous version, multiple connectivity matrices after split
                #list_epoched_data = np.array_split(data,nb_splits,axis = 1)
                
                #print len(list_epoched_data)
                
                #conmat_files = []
                
                #for i,epoched_data in enumerate(list_epoched_data):
                    
                    #conmat_file = compute_and_save_coherency_spectral_connectivity(data=epoched_data,con_method=con_method,sfreq=sfreq,fmin = freq_band[0],fmax = freq_band[1],index = i)
                        
                    #conmat_files.append(conmat_file)
                    
                #return conmat_files
            
############################## testing (can be removed from package) ###########################################

def test_spectral_connectivity():
	
	from params import freq_bands,freq_band_names,con_method
	from mne.connectivity import spectral_connectivity

	subj_path = os.path.join(main_path ,'balai')

	print subj_path

	fif_files = [f for f in os.listdir(subj_path) if f.endswith("fif")]

	print fif_files

	for fif_f in fif_files:

		basename = os.path.splitext(fif_f)[0]

		raw = RawFIF(os.path.join(subj_path,fif_f),preload = True)

		print raw

		print len(raw.ch_names)

		sfreq = raw.info['sfreq']

		select_sensors, = np.where(np.array([ch_name[0] == 'M' for ch_name in raw.ch_names],dtype = 'bool') == True)

		### save electrode locations
		sens_loc = [raw.info['chs'][i]['loc'][:3] for i in select_sensors]
		sens_loc = np.array(sens_loc)

		loc_filename = os.path.join(subj_path,basename +"_correct_channel_coords.txt")
		np.savetxt(loc_filename,sens_loc , fmt = '%s')

		print sens_loc

		### save electrode names

		sens_names = np.array([raw.ch_names[pos] for pos in select_sensors],dtype = "str")
		names_filename = os.path.join(subj_path,basename +"_correct_channel_names.txt")
		np.savetxt(names_filename,sens_names , fmt = '%s')

		#start, stop = raw.time_as_index([0, 100])

		data,times = raw[select_sensors,:]
		print np.max(data,axis = 0)

		for i,freq_band in enumerate(freq_band_names):
			con_matrix, freqs, times, n_epochs, n_tapers = spectral_connectivity(data.reshape(1,data.shape[0],data.shape[1]), method=con_method, mode='multitaper', sfreq=sfreq, fmin= freq_bands[i][0], fmax=freq_bands[i][1], faverage=True, tmin=None,    mt_adaptive=False, n_jobs=1)

			#print con

			con_matrix = np.array(con_matrix[:,:,0])
			print con_matrix.shape
			print np.min(con_matrix),np.max(con_matrix)

			#0/0

			#print data_filtered.shape

			#print data-data
			#print np.max(data-data_filtered,axis = 0)
			#0/0
			np_filename = os.path.join(subj_path,basename+ "_" + con_method +"_" + freq_band +".npy")

			np.save(np_filename,con_matrix)

			#0/0
			
    
if __name__ == '__main__':
    
	test_spectral_connectivity()
