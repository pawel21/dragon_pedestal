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


f = File("Run028.000.fits.fz")
i = 0
event = next(f.Events)
#numer of Events: 7879
#while next(f.Events):
#    event = next(f.Events)
#    i = i + 1
#    print(i)

N = 1000
hiGain_Ch1 = np.zeros(40*N)
for i in range(0, N):
    event = next(f.Events)
    hiGain, lowGain = translate_fits(event)
    hiGain_Ch1[i*40:(i+1)*40] = hiGain[0:40]
    print("i = {} = {}".format(i, hiGain[0:40]))

plt.hist(hiGain_Ch1, bins=50)
plt.show()