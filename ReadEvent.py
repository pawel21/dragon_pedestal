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
size4drs = 4*1024

hi_gain_samples = np.zeros((N, n_channels, RoI))
low_gain_samples = np.zeros((N, n_channels, RoI))
channel_pedestal = np.zeros((n_channels, 2, size4drs))
number_pedestal_event = np.zeros((n_channels, 2, size4drs))

for i in range(0, N):
    event = next(f.Events)
    hiGain, lowGain = translate_fits(event)
    fc = read_first_capacitor(event)
    for j in range(0, n_channels):
        hi_gain_samples[i, j, 0:40] = hiGain[j*40:(j+1)*40]
        low_gain_samples[i, j, 0:40] = lowGain[j*40:(j+1)*40]

    #mean_pedestal for high gain
    for j in range(0, n_channels):
        for k in range(2, RoI-2):
            position = int((k+fc[0][j])%size4drs)
            channel_pedestal[j, 0, position] += hi_gain_samples[i, j, k]
            number_pedestal_event[j, 0, position] += 1


mean_pedestal_high_gain = channel_pedestal[:, 0, :]/number_pedestal_event[:, 0, :]


f2 = File("Run019_Periodic1kHz.fits.fz")
N2 = 400
hi_gain_samples2 = np.zeros((N2, n_channels, RoI))
first_capacitor = np.zeros((N2, 2, 8))
for i in range(0, N2):
    event = next(f.Events)
    hiGain, lowGain = translate_fits(event)
    fc = read_first_capacitor(event)
    for j in range(0, n_channels):
        hi_gain_samples2[i, j, 0:40] = hiGain[j*40:(j+1)*40]
        first_capacitor[i, 0, :] = fc[0, :]


plt.figure()
for i in range(0, 8):
    plt.subplot(2, 4, i+1)
    plt.hist(hi_gain_samples2[:, i, 2:38].flatten(), bins=50)
    std = np.std(hi_gain_samples2[:, i, 2:38].flatten())
    plt.title("Hi Gain Channel {}, $\sigma$ = {:.2f}.".format((i+1), std))


for i in range(0, N2):
    for j in range(0, n_channels):
        for k in range(0, RoI):
            position = int((k + first_capacitor[i, 0, j]) % size4drs)
            hi_gain_samples2[i, j, k] -= mean_pedestal_high_gain[j, position]

plt.figure()
for i in range(0, 8):
    plt.subplot(2, 4, i+1)
    plt.hist(hi_gain_samples2[:, i, 2:38].flatten(), bins=50)
    std = np.std(hi_gain_samples2[:, i, 2:38].flatten())
    plt.title("Hi Gain Channel {}, $\sigma$ = {:.2f}.".format((i + 1), std))

plt.show()