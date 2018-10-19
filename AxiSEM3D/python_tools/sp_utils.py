import numpy as np


def cosine_taper(num_steps):
    
    taper = np.ones(num_steps, dtype = np.float32)

    cut = int(5 * num_steps / 100.)
    cos_part = np.cos(3*np.pi/2. + np.pi/2. * ( np.asarray(range(cut), dtype = np.float32) /(cut-1)))
    taper[0:cut] = cos_part
    taper[-cut:] = np.flip(cos_part,0)

    return taper

def time_to_freq(time):

    freq = 2 * np.pi * np.asarray(range(int(len(time)/2)+1),dtype = np.float32) / (time[-1] - time[0])
    return freq

def log_gabor_filter(freq, center, sigma):

    filt = np.zeros(len(freq), dtype = np.float32)
    filt[1:] = np.exp( -np.power(np.log(np.asarray(freq[1:])/center), 2.)  / (2. * np.power(np.log(sigma),2.)))
    return filt


def pad(trace, num_pad_start, num_pad_end):
    
    return np.pad(trace, (num_pad_start,num_pad_end), 'constant', constant_values=(0,0))
