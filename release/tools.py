import ctypes
import matplotlib.pyplot as plt
import numpy as np


lib = ctypes.cdll.LoadLibrary("./translatezFits.so")
fun = lib.translatezFits
fun.restype = ctypes.c_int
fun.argtypes = [ctypes.POINTER(ctypes.c_ushort),
                ctypes.POINTER(ctypes.c_ushort),
                ctypes.POINTER(ctypes.c_ushort),
                ctypes.POINTER(ctypes.c_ushort),
                ctypes.POINTER(ctypes.c_ushort),
                ctypes.POINTER(ctypes.c_ushort)]


def translate_fits(event):
    eventWfSrcHG = event.hiGain.waveforms.samples
    eventWfSrcLG = event.loGain.waveforms.samples
    eventDRSSrcHG = event.drsTagsHiGain.data
    eventDRSSrcLG = event.drsTagsLoGain.data
    eventWfDestHG = np.zeros((6080,), dtype=np.uint16)
    eventWfDestLG = np.zeros((6080,), dtype=np.uint16)
    nbModules = 19
    c_ushort_p = ctypes.POINTER(ctypes.c_ushort)
    eventWfSrcHG_p = eventWfSrcHG.ctypes.data_as(c_ushort_p)
    eventWfSrcLG_p = eventWfSrcLG.ctypes.data_as(c_ushort_p)
    eventDRSSrcHG_p = eventDRSSrcHG.ctypes.data_as(c_ushort_p)
    eventDRSSrcLG_p = eventDRSSrcLG.ctypes.data_as(c_ushort_p)
    eventWfDestHG_p = eventWfDestHG.ctypes.data_as(c_ushort_p)
    eventWfDestLG_p = eventWfDestLG.ctypes.data_as(c_ushort_p)
    fun(eventWfSrcHG_p, eventWfSrcLG_p, eventDRSSrcHG_p, eventDRSSrcLG_p, eventWfDestHG_p, eventWfDestLG_p, nbModules)
    return eventWfDestHG, eventWfDestLG


def plot_hist(event_data, gain):
    position = [(0,0), (0, 1), (0,2), (0,3),
                (1,0), (1, 1), (1,2), (1,3)]
    fig, ax = plt.subplots(2, 4)
    for i in range(0, 7):
        ax[position[i]].hist(event_data[:, i, 2:38].ravel(), bins=50, facecolor='green', alpha=0.75, histtype='stepfilled')
        ax[position[i]].set_xlabel("signal [counts]")
        ax[position[i]].set_ylabel("number of events")
        sigma = np.std(event_data[:, i, 2:38].ravel())
        mu = np.mean(event_data[:, i, 2:38].ravel())
        ax[position[i]].set_title("{} Channel {}".format(gain, i + 1))
        textstr = '$\mu=%.2f$\n $\sigma=%.2f$' % (mu, sigma)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax[position[i]].text(0.05, 0.95, textstr, transform=ax[position[i]].transAxes, fontsize=14, verticalalignment='top', bbox=props)
    ax[position[7]].axis('off')
    plt.subplots_adjust(top=0.925, bottom=0.08, left=0.080, right=0.95, hspace=0.305, wspace=0.440)