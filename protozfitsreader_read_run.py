from protozfits.simple import File
from enum import Enum
import matplotlib.pyplot as plt


f = File("Run019_Periodic1kHz.fits.fz")

for row in f.RunHeader:
    pass

row._asdict()


print(f.Events.header)
#32678 / Offset for uint16

event = next(f.Events)
event._asdict()
print(event._fields)

#higain = event.hiGain.waveforms.samples

print(event.cameraCounters)
print(event.cameraCounters.counters)

#lo_gain_samples = toNumPyArray(
#    event.loGain.waveforms.samples
#).reshape(-1, event.loGain.waveforms.num_samples)

t = []
time = []
oldtime = []

for i in range(0, 100):
    event = next(f.Events)
    word1 = event.cameraCounters.counters[7]
    word2 = event.cameraCounters.counters[8]
    word3 = event.cameraCounters.counters[9]
    word4 = event.cameraCounters.counters[10]
    newtime = (word1 + word2*2**16 + word3*2**32 + word4*2**48)/(133.e6)
    tt = (word1 + word2*2**16 + word3*2**32 + word4*2**48)/(133.e6)
    time.append(tt)
    if i > 0:
        t.append(newtime-oldtime)
    oldtime = newtime

plt.plot(time)
plt.show()

print(event.cameraCounters.counters[:11])

eventWfSrcHG = event.hiGain.waveforms.samples
eventWfSrcLG = event.loGain.waveforms.samples
eventDRSSrcHG=  event.drsTagsHiGain.data
eventDRSSrcLG = event.drsTagsLoGain.data
# descination
eventWfDestHG = []
eventWfDestLG = []
nbModules = 19


def translateFits():
    i, j, k = 0, 0, 0
