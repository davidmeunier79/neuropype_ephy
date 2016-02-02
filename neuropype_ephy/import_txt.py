# -*- coding: utf-8 -*-


def split_txt(sample_size,txt_file,sep_label_name):

    import os

    import numpy as np
    import pandas as pd

    df = pd.read_table(txt_file,sep = ";",decimal = ",", header = None, index_col = 0)


    ## electrode names:
    np_indexes = df.index.values

    list_indexes = np_indexes.tolist()

    print list_indexes

    if sep_label_name != "" :
        keep = np.array([len(index.split(sep_label_name)) == 2 for index in list_indexes],dtype = "int")
        
    else:
        keep = np.ones(shape = np_indexes.shape)
                        
    print keep

    print np_indexes[keep == 1]


    elec_names_file = os.path.abspath("correct_channel_names.txt")

    np.savetxt(elec_names_file,np_indexes[keep == 1],fmt = "%s")

    ## splitting data_path
    print df.shape

    if df.shape[1] % sample_size != 0:

        print "Error, sample_size is not a multiple of ts shape"

        print sample_size
        print df
        
        0/0
        
        return 

    nb_epochs = df.shape[1] / sample_size

    print nb_epochs

    splitted_ts = np.split(df.values,nb_epochs,axis = 1)

    print len(splitted_ts)

    print splitted_ts[0]

    np_splitted_ts = np.array(splitted_ts,dtype = 'float')

    print np_splitted_ts.shape

    splitted_ts_file = os.path.abspath("splitted_ts.npy")

    np.save(splitted_ts_file,np_splitted_ts)

    return splitted_ts_file,elec_names_file


def split_txt_multiwin(sample_size,txt_file,n_windows,sep_label_name):

    #, mode = 'multiplexed'
    import os

    import numpy as np
    import pandas as pd

    #if mode == 'multiplexed':
        #df = pd.read_table(txt_file,sep = ";",decimal = ",", header = 0, index_col = None)
    #else:
        #df = pd.read_table(txt_file,sep = ";",decimal = ",", header = None, index_col = 0)
        
    df = pd.read_table(txt_file,sep = ";",decimal = ",", header = None, index_col = 0)
    print df.values.shape
    
    ## electrode names:
    np_indexes = df.index.values

    list_indexes = np_indexes.tolist()

    print list_indexes

    if sep_label_name != "" :
        keep = np.array([len(index.split(sep_label_name)) == 2 for index in list_indexes],dtype = "int")        
    else:
        keep = np.ones(shape = np_indexes.shape)
                        
    print keep

    print np_indexes[keep == 1]


    elec_names_file = os.path.abspath("correct_channel_names.txt")

    np.savetxt(elec_names_file,np_indexes[keep == 1],fmt = "%s")

    ## splitting data_path
    print df.shape

    if df.shape[1] % sample_size != 0:

        print "Error, sample_size is not a multiple of ts shape"

        print sample_size
        print df
        
        0/0
        
        return 
    
    nb_epochs = df.shape[1] / sample_size

    print nb_epochs

    splitted_ts = np.split(df.values,nb_epochs,axis = 1)

    print len(splitted_ts)

    print splitted_ts[0]

    np_splitted_ts = np.array(splitted_ts,dtype = 'float')

    print np_splitted_ts.shape

    print len(n_windows)
    
    win_splitted_ts = np.array([np_splitted_ts[:,:,n_win[0]:n_win[1]] for n_win in n_windows])
    
    #win_splitted_ts = []
    
    #for n_win in n_windows:
        
        #win_splitted_ts.append(np_splitted_ts[:,:,n_win[0]:n_win[1]])
        
    #np_win_splitted_ts = np.array(win_splitted_ts)
    
    
    print win_splitted_ts.shape
    
    win_splitted_ts_file = os.path.abspath("win_splitted_ts.npy")

    np.save(win_splitted_ts_file,win_splitted_ts)

    return win_splitted_ts_file,elec_names_file
