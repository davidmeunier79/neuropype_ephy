"""Definitions of multiscale entropy wrappers"""

def get_mse_multiple_sensors(ts_file, m, r):
    """Compute multiscale entropy across sensors for epoched data
    ep_data should be in {n_trials x n_sensors x n_times} format

    """
    import numpy as np
    import os
    from neuropype_ephy.mse import get_mse_curve_across_trials
    ep_data = np.load(ts_file)
    mean_mses = []
    std_mses = []
    for isen in range(ep_data.shape[1]):
        single_sensor_data = ep_data[:,isen,:].T
        mean_mse, std_mse = get_mse_curve_across_trials(single_sensor_data, m, r)
        mean_mses.append(mean_mse)
        std_mses.append(std_mse)
    mean_mses = np.array(mean_mses)
    std_mses = np.array(std_mses)
    savename = os.path.abspath('mse.npz')
    np.savez(savename, mean_mse=mean_mses, std_mse=std_mses)
    return savename
        


def get_mse_curve_across_trials(single_sensor_data, m, r):
    """Assumption: single_sensor_data is in (time trial) format
    mse is column vector containing mean entropy at different scales
    """

    import numpy as np
    import os
    import glob

    dir_path = os.path.dirname(os.path.realpath(__file__))
    out_file = open('data_filelist.txt', 'w')
    for itrial in range(single_sensor_data.shape[1]):
        data = single_sensor_data[:, itrial]
        fname = 'data_trial' + str(itrial + 1) + '.txt'
        np.savetxt(fname, data)
        out_file.write('{} \n'.format(fname))
    out_file.close()

    cmd = dir_path + '/mse -n 40 -a 1 -m ' + str(m) + ' -M ' +  str(m) +\
          ' -r ' + str(r) + ' -R ' + str(r) + ' -F data_filelist.txt > data_filelist.mse'
    os.system(cmd)
    res = np.loadtxt('data_filelist.mse', skiprows=7)
    flist = glob.glob('data_trial*.txt')
    for fname in flist:
        os.remove(fname)
    os.remove('data_filelist.txt')
    os.remove('data_filelist.mse')

    mean_mse = res[:, 1]
    std_mse = res[:, 2]

    return mean_mse, std_mse



if __name__ == '__main__':
    test_path = '/media/dmalt/SSD500/aut_gamma/aut_gamma_pipeline/\
_keys_K0001__K0001ec-epo-fif/ep2ts/ts_epochs.npy'
    import numpy as np
    # data1 = data[:,0,:] 
    m = 2 
    r = 0.2

    # mean_mse, std_mse = get_mse_curve_across_trials(data1.T, m, r)
    get_mse_multiple_sensors(test_path, m, r)
    # print(mean_mse)
    # print(std_mse)

