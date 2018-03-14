import matplotlib.pyplot as plt
import numpy as np
from protozfits.simple import File
from ctapipe.core import Container, Field
from enum import Enum


class PartialEventCounters:
    def __init__(self, event_counters):
        self.event_counters = event_counters
        self.ppsCounter = np.uint16(0)
        self.tenMhzCounter = np.uint32(0)
        self.absEventCounter = np.uint32(0)
        self.triggerCounter = np.uint32(0)
        self.localClockCounter = np.uint64(0)

    def read_event_counters(self):
        self.ppsCounter = self.event_counters[0]
        self.tenMhzCounter = self.event_counters[1] + self.event_counters[2]*2**16
        self.absEventCounter = self.event_counters[3] + self.event_counters[4]*2**16
        self.triggerCounter = self.event_counters[5] + self.event_counters[6]*2**16
        self.localClockCounter = self.event_counters[7] + self.event_counters[8]*2**16 +   \
                                 self.event_counters[9]*2**32 + self.event_counters[10]*2**48


f = File("Run019_Periodic1kHz.fits.fz")
event = next(f.Events)

list_of_events = []
len_event = int(len(event.cameraCounters.counters)/11.0)
for i in range(0, len_event):
    e = event.cameraCounters.counters[i*11:i*11+11]
    list_of_events.append(PartialEventCounters(e))
    list_of_events[i].read_event_counters()
    print(i, e)
