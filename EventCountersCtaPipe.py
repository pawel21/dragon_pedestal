import numpy as np
from protozfits.simple import File
from ctapipe.core import Container, Field


class EventCountersContainer(Container):
    ppsCounter = Field(np.uint16, "np.uint16")
    tenMhzCounter = Field(np.uint32, "np.uint32")
    absEventCounter = Field(np.uint32, "np.uint32")
    triggerCounter = Field(np.uint32, "np.uint32")
    localClockCounter = Field(np.uint64, "np.unit64")


def create_stream(n_event, event_counters):
    data = EventCountersContainer()
    for i in range(n_event):
        data.ppsCounter = event_counters[0]
        data.tenMhzCounter = event_counters[1] + event_counters[2] * 2 ** 16
        data.absEventCounter = event_counters[3] + event_counters[4] * 2 ** 16
        data.triggerCounter = event_counters[5] + event_counters[6] * 2 ** 16
        data.localClockCounter = event_counters[7] + event_counters[8] * 2 ** 16 + \
                         event_counters[9] * 2 ** 32 + event_counters[10] * 2 ** 48

        yield data


f = File("Run019_Periodic1kHz.fits.fz")
event = next(f.Events)
e = event.cameraCounters.counters[0:11]

for data in create_stream(1, e):
    for key, val in data.items():
        print("{}: {}, type : {}".format(key, val, type(val)))