import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
from protozfits.simple import File
from enum import Enum


#read file
f = File("Run021.1.fits.fz")
event = next(f.Events)

eventWfSrcHG = event.hiGain.waveforms.samples
eventWfSrcLG = event.loGain.waveforms.samples
eventDRSSrcHG=  event.drsTagsHiGain.data
eventDRSSrcLG = event.drsTagsLoGain.data
# descination
eventWfDestHG = np.zeros(eventWfSrcHG.shape, dtype=np.uint16)
eventWfDestLG = np.zeros(eventWfSrcLG.shape, dtype=np.uint16)
nbModules = 10

# load c libary to use c code
lib = ctypes.cdll.LoadLibrary("./translatezFits.so")
fun = lib.translatezFits
fun.restype = ctypes.c_int
fun.argtypes = [ctypes.POINTER(ctypes.c_ushort),
                ctypes.POINTER(ctypes.c_ushort),
                ctypes.POINTER(ctypes.c_ushort),
                ctypes.POINTER(ctypes.c_ushort),
                ctypes.POINTER(ctypes.c_ushort),
                ctypes.POINTER(ctypes.c_ushort)]

c_ushort_p = ctypes.POINTER(ctypes.c_ushort)
eventWfSrcHG_p = eventWfSrcHG.ctypes.data_as(c_ushort_p)
eventWfSrcLG_p = eventWfSrcLG.ctypes.data_as(c_ushort_p)
eventDRSSrcHG_p = eventDRSSrcHG.ctypes.data_as(c_ushort_p)
eventDRSSrcLG_p = eventDRSSrcLG.ctypes.data_as(c_ushort_p)
eventWfDestHG_p = eventWfDestHG.ctypes.data_as(c_ushort_p)
eventWfDestLG_p = eventWfDestLG.ctypes.data_as(c_ushort_p)

fun(eventWfSrcHG_p, eventWfSrcLG_p, eventDRSSrcHG_p, eventDRSSrcLG_p, eventWfDestHG_p, eventWfDestLG_p, nbModules)