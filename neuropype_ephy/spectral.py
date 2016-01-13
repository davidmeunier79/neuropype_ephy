# -*- coding: utf-8 -*-

import numpy as np


################################################### compute spectral connectivity #############################################################################"

def compute_and_save_spectral_connectivity(data,con_method,sfreq,fmin,fmax,index = 0,mode = 'multitaper'):

    import sys,os
    from mne.connectivity import spectral_connectivity

    import numpy as np

    print data.shape

    if len(data.shape) < 3:
        if con_method in ['coh','cohy','imcoh']:
            data = data.reshape(1,data.shape[0],data.shape[1])

        elif con_method in ['pli','plv','ppc' ,'pli','pli2_unbiased' ,'wpli' ,'wpli2_debiased']:
            print "warning, only work with epoched time series"
            sys.exit()
        
    if mode == 'multitaper':
        
        con_matrix, freqs, times, n_epochs, n_tapers  = spectral_connectivity(data, method=con_method, sfreq=sfreq, fmin= fmin, fmax=fmax, faverage=True, tmin=None, mode = 'multitaper',   mt_adaptive=False, n_jobs=1)
        
        con_matrix = np.array(con_matrix[:,:,0])

    elif mode == 'cwt_morlet':
        
        frequencies = np.arange(fmin, fmax, 1)
        n_cycles = frequencies / 7.

        con_matrix, freqs, times, n_epochs, n_tapers  = spectral_connectivity(data, method=con_method, sfreq=sfreq, faverage=True, tmin=None, mode='cwt_morlet',   cwt_frequencies= frequencies, cwt_n_cycles= n_cycles, n_jobs=1)
        
        con_matrix = np.mean(np.array(con_matrix[:,:,0,:]),axis = 2)
    
    else:
        
        print "Error, mode = %s not implemented"%(mode)
        
        return []

    print con_matrix.shape
    print np.min(con_matrix),np.max(con_matrix)

    conmat_file = os.path.abspath("conmat_" + str(index) + "_" + con_method + ".npy")

    np.save(conmat_file,con_matrix)

    return conmat_file

def plot_circular_connectivity(conmat_file, labels_file, is_sensor_space, nb_lines = 200):

    import os
    
    import numpy as np
    
    from nipype.utils.filemanip import split_filename as split_f
    
    from mne.viz import circular_layout, plot_connectivity_circle
    import matplotlib.pyplot as plt
    
    if is_sensor_space:
        label_names = [line.strip() for line in open(labels_file)]
        node_order  = label_names
        node_colors = None
        
    else:
        node_colors = [label.color for label in labels_file]    
        # reorder the labels based on their location in the left hemi
        label_names = [label.name for label in labels_file]

        lh_labels = [name for name in label_names if name.endswith('lh')]

        # Get the y-location of the label
        label_ypos = list()
        for name in lh_labels:
            idx = label_names.index(name)
            ypos = np.mean(labels_file[idx].pos[:, 1])
            label_ypos.append(ypos)

        # Reorder the labels based on their location
        lh_labels = [label for (yp, label) in sorted(zip(label_ypos, lh_labels))]
        
        # For the right hemi
        rh_labels = [label[:-2] + 'rh' for label in lh_labels]

        # Save the plot order 
        node_order = list()
        node_order.extend(lh_labels[::-1])  # reverse the order
        node_order.extend(rh_labels)

    
    path,fname,ext = split_f(conmat_file)    
    print fname
    
    
    #print label_names

    conmat = np.load(conmat_file)
    print conmat.shape
    
    
    # Angles
    node_angles = circular_layout(label_names, node_order, start_pos=90,
                                group_boundaries=[0, len(label_names) / 2])

    # Plot the graph using node colors from the FreeSurfer parcellation. We only
    # show the 300 strongest connections.
    fig,_ = plot_connectivity_circle(conmat, label_names, n_lines=nb_lines,  
                                     node_angles=node_angles, node_colors = node_colors,
                                     fontsize_names = 12, 
                                     title='All-to-All Connectivity' , show = False)
    
    
    #plot_conmat_file = os.path.abspath('circle.png')
    plot_conmat_file = os.path.abspath('circle_' + fname + '.eps')
    fig.savefig(plot_conmat_file, facecolor='black')
    
    
    plt.close(fig)
    #fig1.close()
    del fig
    
    return plot_conmat_file
    
#################################################################################################################################################################"

def spectral_proc(ts_file,sfreq,freq_band,con_method):

    import numpy as np
    #import os

    from neuropype_ephy.spectral import compute_and_save_spectral_connectivity

    data = np.load(ts_file)

    conmat_file = compute_and_save_spectral_connectivity(data = data,con_method = con_method,sfreq=sfreq, fmin= freq_band[0], fmax=freq_band[1])
        
    return conmat_file


def spectral_proc_label(ts_file,sfreq,freq_band,con_method,label,mode):

    import numpy as np
    #import os

    from neuropype_ephy.spectral import compute_and_save_spectral_connectivity

    data = np.load(ts_file)

    conmat_file = compute_and_save_spectral_connectivity(data = data,con_method = con_method,sfreq=sfreq, fmin= freq_band[0], fmax=freq_band[1],index = label,mode = mode)

    return conmat_file


def multiple_spectral_proc(ts_file,sfreq,freq_band,con_method):

    import numpy as np
    import os

    from mne.connectivity import spectral_connectivity

    all_data = np.load(ts_file)

    print all_data.shape
    
    #print sfreq
            
    print freq_band
    
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


def epoched_spectral_proc(ts_file,sfreq,freq_band,freq_band_name,con_method,epoch_window_length):

    import numpy as np

    from neuropype_ephy.spectral import compute_and_save_spectral_connectivity

    data = np.load(ts_file)

    print data.shape
    print sfreq
    print freq_band
    print freq_band_name

    if epoch_window_length == None:
        
        conmat_file = compute_and_save_spectral_connectivity(data=data,con_method=con_method,sfreq=sfreq,fmin = freq_band[0],fmax = freq_band[1])
    else:
            
            print "Shape before splits:"
            print data.shape
            
            print  "sfreq:"
            print sfreq
            
            nb_splits = data.shape[1] // (epoch_window_length * sfreq)
            
            print "nb_splits:"
            print nb_splits
            
            reste = data.shape[1] % int(epoch_window_length * sfreq)
            
            print "reste:"
            print reste
            
            if reste != 0:
                data = data[:,:-reste]
            
            print "shape after reste:"
            print data.shape
            
            print "epoching data with {}s by window, resulting in {} epochs".format(epoch_window_length,nb_splits)
            
            
            
            list_epoched_data = np.array_split(data,nb_splits,axis = 1)
            
            for epo in list_epoched_data:
                print epo.shape
            
            #print "Shape after splits:"
            #print epoched_data.shape

            epoched_data = np.array(list_epoched_data)
            
            conmat_file = compute_and_save_spectral_connectivity(data=epoched_data, con_method=con_method, sfreq=sfreq, fmin= freq_band[0], fmax=freq_band[1])

            return conmat_file
        
        
        
########################################################### plot spectral connectivity #################################################################

def plot_circular_connectivity(conmat_file,labels_file,nb_lines, vmin = None, vmax = None):

    import os
    
    import numpy as np
    
    from nipype.utils.filemanip import split_filename as split_f
    
    from mne.viz import circular_layout, plot_connectivity_circle
    import matplotlib.pyplot as plt
    
    label_names= [line.strip() for line in open(labels_file)]
    
    path,fname,ext = split_f(conmat_file)
    
    print fname
    
    
    #print label_names
    conmat = np.load(conmat_file)
    #print conmat.shape
    
    # Angles
    node_angles = circular_layout(label_names, node_order = label_names, start_pos=90,
                                group_boundaries=[0, len(label_names) / 2])

    # Plot the graph using node colors from the FreeSurfer parcellation. We only
    # show the 300 strongest connections.
    fig,_ = plot_connectivity_circle(conmat, label_names, n_lines=nb_lines,  node_angles=node_angles, fontsize_names = 12, title='All-to-All Connectivity' , show = False, vmin = vmin, vmax = vmax)
    
    
    #plot_conmat_file = os.path.abspath('circle.png')
    plot_conmat_file = os.path.abspath('circle_' + fname + '.eps')
    fig.savefig(plot_conmat_file, facecolor='black')
    
    
    plt.close(fig)
    #fig1.close()
    del fig
    
    return plot_conmat_file

     
def filter_adj_plot_mat(conmat_file,labels_file,sep_label_name,k_neigh):

    import numpy as np
    import os
    
    from itertools import combinations
    
    labels = [line.strip().split(sep_label_name) for line in open(labels_file)]
    
    print labels
    
    triu_indices = np.triu_indices(len(labels),1)
    
    print triu_indices
            
    adj_mat = np.zeros(shape = (len(labels),len(labels)),dtype = bool)
                   
    for i in range(k_neigh):               
                    
        adj_plots = [(a[0] == b[0]) and ((int(a[1]) + i + 1 )== int(b[1])) for a,b in combinations(labels,2)]
        
        print len(adj_plots)
        
        adj_mat[triu_indices] =  adj_mat[triu_indices] + adj_plots
        
        print np.sum(adj_mat == True)
        
    print adj_mat
    
    ### loading ad filtering conmat_file
    conmat = np.load(conmat_file)
    
    print conmat
    
    assert conmat.shape[0] == len(labels), "warning, wrong dimensions between labels and conmat"
    
    filtered_conmat = np.transpose(conmat).copy()
    
    #print np.transpose(adj_mat) == True
    
    x,y = np.where(adj_mat == True)
    
    filtered_conmat[x,y] = 0.0
    
    print filtered_conmat
    
    filtered_conmat_file = os.path.abspath("filtered_conmat.npy")
    
    np.save(filtered_conmat_file,np.transpose(filtered_conmat))
    
    return filtered_conmat_file
  
            
            
            
            
            
            
            
############################## testing (can be removed from package) ###########################################

def test_spectral_connectivity(main_path = "/mnt/Data/Projet-Karim", con_method = 'coh',  freq_bands = [[15.,40.]], freq_band_names = ['beta']):

    import os 
    from mne.connectivity import spectral_connectivity

    from mne.io import RawFIF
    
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
