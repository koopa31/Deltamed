import numpy as np
import matcompat
import os

from io import StringIO
# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass


def import_deltamed(filename):

    """
    import_deltamed allows to load a Deltamed signal loading the metadata from the .txt file and loading the signal
    from a binary or an ASCII file (loading the binary is much more efficient from a computational time point of view)

    :param filename: name of the txt file
    :return: [out_header, out_data, message_string] list of the header (metadata), the signal, and a string containing
    warning messages.
    """
    message_string = []
    out_header = np.array([])
    out_data = np.array([])
    event_pos = []
    event_code = []
    message_string.append(np.array(np.hstack(('Loading : ', filename))))
    # filename_txt
    p = os.path.split(filename)[0]
    n = os.path.split(filename)[1].split('.')[0]
    e = '.' + os.path.split(filename)[1].split('.')[1]
    filename_txt = os.path.join(p, n) + e
    message_string.append(np.array(np.hstack(('Loading TXT header : ', filename_txt))))
    # open text header
    txt_file = open(filename_txt, encoding='iso8859_14')
    tp = txt_file.readline().rstrip()
    while type(tp) == str:
        if len(tp) >= 8.:
            if tp[0:8] in 'Sampling':
                sampling_rate = float(tp[9:len(tp)])
                message_string[int((0+1.))-1] = np.array(np.hstack(('Sampling rate : ', str(sampling_rate))))

            if tp[0:8] in 'Channels':
                st = tp[9:len(tp)]

            if tp[0:8] in 'Gainx100':
                st2 = tp[10:len(tp)]
        
        if len(tp) >= 7.:
            if tp[0:7] in '[EVENT]':
                break

        tp = txt_file.readline().rstrip()
        
    index = 1.
    while type(tp) is str and len(tp) != 0:
        tp = txt_file.readline().rstrip()
        if type(tp) is str and len(tp) != 0:
            idx = tp.find(',')
            event_pos.append(int(tp[0:idx]))
            event_code.append(tp[int(idx+1):])
            index = index + 1

    channel_labels = np.loadtxt(StringIO(st), delimiter=',', dtype=str)
    channel_gain = np.loadtxt(StringIO(st2), delimiter=',', dtype=float)
    xstart = 0.
    xstep = 1./sampling_rate
    event_lat = (np.double(event_pos)-1.)*xstep
    # generate header
    out_header = {}
    out_header['filetype'] = 'time amplitude'
    out_header['name'] = filename_txt
    out_header['tags'] = ''
    # header.history
    out_header['history'] = np.array([])
    # initialize header.datasize (will be updated later)
    out_header['datasize'] = np.array(np.hstack((0., 0., 0., 0., 0., 0.)))
    out_header['xstart'] = 0.
    out_header['ystart'] = 0.
    out_header['zstart'] = 0.
    out_header['xstep'] = xstep
    out_header['ystep'] = 1.
    out_headerstep = 1
    # chanlocs
    out_header['chanlocs'] = {}
    out_header['chanlocs']['labels'] = []
    out_header['chanlocs']['topo_enabled'] = []
    for i in np.arange(1, (len(channel_labels))+1):
        out_header['chanlocs']['labels'].append(channel_labels[int(i)-1])
        out_header['chanlocs']['topo_enabled'].append(0)
        
    message_string[int((0+1))-1] = np.array(np.hstack(('Number of channels : ', str(len((out_header['chanlocs']))))))
    # filename_bin
    filename_bin = os.path.join(p, n) + '.bin'
    # load BIN data if it exists, else will try to load ASC data
    # note that Deltamed software can export as either BIN or ASCII
    if os.path.exists(filename_bin):
        tp1 = np.fromfile(filename_bin, 'int16')
        epoch_size = int(len(tp1)/len(out_header['chanlocs']['labels']))
        message_string[0] = np.array(np.hstack(('Epoch size : ', str(epoch_size))))
        tp2 = np.reshape(tp1, (len((out_header['chanlocs']['labels'])), epoch_size))
        out_data = np.zeros(shape=(1, tp2.shape[0], 1, 1, 1, tp2.shape[1]))
        for chanpos in np.arange(0, tp2.shape[0]):
            out_data[0, int(chanpos)-1, 0, 0, 0, :] = np.dot(tp2[int(chanpos)-1, :], channel_gain[int(chanpos)-1])/1000.
            
    else:
        # filename_asc
        filename_asc = os.path.join(p, n) + ".asc"
        # load ASC data it it exists
        if os.path.exists(filename_asc):
            message_string[int((0+1))-1] = 'Loading the ASCII data. *** This can take a while ***'
            message_string[int((0+1))-1] = 'For faster operations, try to export data as BIN'
            tp = np.loadtxt(filename_asc, encoding='ascii')
            out_data = np.zeros(shape=(1, len(np.arange(1., (matcompat.size(tp, 2.))+1)), 1, 1, 1,
                                       len(np.dot(tp[:, 0], channel_gain[0])/1000.)))
            for i in np.arange(1., (matcompat.size(tp, 2.))+1):
                out_data[0, int(i) - 1, 0, 0, 0, :] = np.dot(tp[:, int(i)-1], channel_gain[int(i)-1])/1000.
                
        else:
            message_string[0] = 'Error : no data file found!'
            return []

    # update header.datasize
    out_header['datasize'] = matcompat.size(out_data)
    message_string[0] = np.array(np.hstack(('Number of epochs : ', str((out_header['datasize'][0])))))
    return [out_header, out_data, message_string]
