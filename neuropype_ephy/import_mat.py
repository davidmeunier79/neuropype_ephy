# -*- coding: utf-8 -*-

import os
import mne

from mne import filter
#from mne.io import RawFIF

#mne.set_log_level('WARNING')

import numpy as np

#from params import main_path

##################### nodes (Function)

#def preprocess_mat_to_conmat(mat_file,):
def import_mat_to_conmat(mat_file,orig_channel_names_file,orig_channel_coords_file):

    import os
    import numpy as np

    import mne

    from mne.io import RawArray	
    from nipype.utils.filemanip import split_filename as split_f

    from scipy.io import loadmat

    subj_path,basename,ext = split_f(mat_file)

    mat = loadmat(mat_file)

    #field_name = basename.split('_')[0]
    #field_name = basename.split('_')[1]
    #print field_name
    print mat["mat_x"].shape


    raw_data = np.array(mat["mat_x"],dtype = "f")
    print raw_data
    print raw_data.shape

    conmat_file = os.path.abspath(basename +".npy")

    np.save(conmat_file,raw_data)

    correct_channel_coords = np.loadtxt(orig_channel_coords_file)
    
    print correct_channel_coords
    
    
    correct_channel_names = np.loadtxt(orig_channel_coords_file)
    
    print correct_channel_names
    
    ### save channel coords
    channel_coords_file = os.path.abspath("correct_channel_coords.txt")
    np.savetxt(channel_coords_file ,correct_channel_coords , fmt = '%s')

    ### save channel names
    channel_names_file = os.path.abspath("correct_channel_names.txt")
    np.savetxt(channel_names_file,correct_channel_names , fmt = '%s')

    return conmat_file,channel_coords_file,channel_names_file

def import_tsmat_to_ts(tsmat_file,data_field_name = 'F', good_channels_field_name = 'ChannelFlag', hdf5_mat = False):
    #,orig_channel_names_file,orig_channel_coords_file):

    import os
    import numpy as np

    import mne

    from mne.io import RawArray	
    from nipype.utils.filemanip import split_filename as split_f

    print tsmat_file
    subj_path,basename,ext = split_f(tsmat_file)

    if hdf5_mat:
        import h5py
        f = h5py.File(tsmat_file, 'r')
        raw_data=np.asarray(list(f[data_field_name])).T
    else:
        from scipy.io import loadmat
        mat = loadmat(tsmat_file)
        raw_data = np.array(mat[data_field_name],dtype = "f")

    print raw_data.shape

    #print [key for key in mat.keys()]
    
    #if good_channels_field_name != None:
        
        #print "Using good channels to sort channels"
        #good_channels = np.array(mat[good_channels_field_name])
        #print good_channels.shape
        
        
        #good_channels = good_channels.reshape(good_channels.shape[0])
        #print good_channels.shape
        
        
        #good_data = raw_data[good_channels == 1,:]
        
        #print good_data.shape
        
    #else:
        #print "No channel sorting" 
        #good_data = raw_data
        
	good_data = raw_data
	
    #### save data
    print good_data.shape
    
    ts_file = os.path.abspath("tsmat.npy")
    np.save(ts_file,good_data)
    return ts_file


def import_amplmat_to_ts(tsmat_file):
    #,orig_channel_names_file,orig_channel_coords_file):

    import os
    import numpy as np

    import mne

    from mne.io import RawArray	
    from nipype.utils.filemanip import split_filename as split_f

    from scipy.io import loadmat

    print tsmat_file
        
    subj_path,basename,ext = split_f(tsmat_file)

    mat = loadmat(tsmat_file)

    #field_name = basename.split('_')[0]
    #field_name = basename.split('_')[1]
    #print field_name

    raw_data = np.array(mat['F'],dtype = "f")
    print raw_data.shape

    good_channels = np.array(mat['ChannelFlag'])

    good_channels = good_channels.reshape(good_channels.shape[0])

    print "Good channels:"
    print good_channels.shape

    good_data = raw_data[good_channels == 1,:]

    print good_data.shape

    #### save data 
    ts_file = os.path.abspath("amplmat.npy")

    np.save(ts_file,good_data)
    #np.save(ts_file,raw_data)

    return ts_file




def import_mat_to_ts(mat_file,orig_channel_names_file,orig_channel_coords_file):

    import os
    import numpy as np

    import mne

    from mne.io import RawArray	
    from nipype.utils.filemanip import split_filename as split_f

    from scipy.io import loadmat

    subj_path,basename,ext = split_f(mat_file)

    mat = loadmat(mat_file)

    field_name = basename.split('_')[0]
#    field_name = basename.split('_')[1]
    print '************************************ \n'
    print field_name
    print mat[field_name].shape


    raw_data = np.array(mat[field_name],dtype = "f")
    print raw_data.shape

    #### load electrode names
    elec_names = [line.strip() for line in open(orig_channel_names_file)]

    print elec_names
    #0/0

    #select_sensors = np.arange(21)

    select_sensors, = np.where(np.array([elec.startswith( 'EMG') or elec.startswith('EOG') for elec in elec_names],dtype = 'bool') == False)
    print select_sensors

    ### save electrode locations
    elec_loc = np.loadtxt(orig_channel_coords_file)
    print elec_loc

    correct_elec_loc = np.roll(elec_loc[select_sensors,:],shift = 2,axis = 1)

    print correct_elec_loc[0,:]
    print correct_elec_loc[7,:]
        
        
    channel_coords_file = os.path.abspath("correct_channel_coords.txt")
    np.savetxt(channel_coords_file ,correct_elec_loc , fmt = '%s')

    ### save electrode names
    correct_elec_names = np.array([elec_names[pos] for pos in select_sensors],dtype = "str")

    channel_names_file = os.path.abspath("correct_channel_names.txt")
    np.savetxt(channel_names_file,correct_elec_names , fmt = '%s')

    ### save data (reorganise dimensions)
    new_data = raw_data[select_sensors,:].swapaxes(0,2).swapaxes(1,2)
    print new_data.shape


    ts_file = os.path.abspath(basename +".npy")

    np.save(ts_file,new_data)

    return ts_file,channel_coords_file,channel_names_file

def concat_ts(all_ts_files):
    
    import numpy as np
    import os
    
    print all_ts_files

    for i,ts_file in enumerate(all_ts_files):
    
        ## loading ROI coordinates
        ts = np.load(ts_file)
        
        #print "all_ts: " 
        print ts.shape
        
        if i == 0:
            concat_ts = ts.copy()
            #print concat_ts.shape
        else:
            concat_ts = np.concatenate((concat_ts,ts),axis = 0)
            #print concat_ts.shape

    print concat_ts.shape

    ### saving time series
    concat_ts_file = os.path.abspath("concat_ts.npy")
    np.save(concat_ts_file,concat_ts)
    
        
    return concat_ts_file

if __name__ == '__main__':
    
    test_loadmat()

# $MNE_ROOT/bin/mne_ctf2fiff --ds balai_speech-oscill_20130716_01.ds/

