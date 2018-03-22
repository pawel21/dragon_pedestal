import ctypes
import matplotlib.pyplot as plt
import numpy as np
from protozfits.simple import File
from enum import Enum


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


#only for modules one
def read_first_capacitor(event):
    first_capacitor = event.firstCapacitorIds.data[0:8]
    high = 0; low = 1
    fc = np.zeros((2, 8))

    for i in [0, 2, 4, 6]:
        fc[high][i] = first_capacitor[i]
        fc[high][i+1] = first_capacitor[i]

    for j in [1, 3, 5, 7]:
        fc[low][j-1] = first_capacitor[j]
        fc[low][j] = first_capacitor[j]

    return fc

f = File("Run021.1.fits.fz")
N = 15000
n_channels = 8
n_modules = 19
RoI = 40
hi_gain_samples = np.zeros((N, n_channels, RoI, n_modules))
low_gain_samples = np.zeros((N, n_channels, RoI, n_modules))
hiGain_Ch1 = np.zeros(40*N)

for i in range(0, N):
    event = next(f.Events)
    hiGain, lowGain = translate_fits(event)
    fc = read_first_capacitor(event)
    for j in range(0, n_modules):
        for k in range(0, n_channels):
            hi_gain_samples[i, k, 0:40, j] = hiGain[(j+k)*40:(j+k+1)*40]
            low_gain_samples[i, k, 0:40, j] = lowGain[(j + k) * 40:(j + k + 1) * 40]


plt.figure()
for i in range(0, 8):
    plt.subplot(2, 4, i+1)
    plt.hist(hi_gain_samples[:, i, 2:38, 1].flatten(), bins=50)
    plt.title("Hi Gain Channel {}".format(i+1))

plt.figure()
for i in range(0, 8):
    plt.subplot(2, 4, i+1)
    plt.hist(low_gain_samples[:, i, 2:38, 1].flatten(), bins=50)
    plt.title("Low Gain Channel {}".format(i+1))

plt.show()