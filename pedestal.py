import numpy as np

from event import Event


class PedestalSimple:
    nch = 8
    roisize = 40
    size4drs = 4 * 1024

    def __init__(self):
        self.meanped = np.zeros((self.nch, 2, self.size4drs), dtype=np.float64)
        self.rmsped = np.zeros((self.nch, 2, self.size4drs), dtype=np.float64)
        self.numped = np.zeros((self.nch, 2, self.size4drs), dtype=np.int32)

    def fill_ped_event(self, evt):
        for i in range(0, self.nch):
            for j in range(0, 2):
                fc = int(evt.firstcap[i,j])
                for k in range(0, evt.roisize-2):
                    posabs = int((k+fc)%self.size4drs)
                    val = evt.samples[i, j, k]
                    self.meanped[i, j, posabs] += val
                    self.rmsped[i, j, posabs] += val*val
                    self.numped[i, j, posabs] += 1

    def finalize_ped(self):
        for i in range(0, self.nch):
            for j in range(0, 2):
                for k in range(0, self.size4drs):
                    if self.numped[i, j, k]:
                        self.meanped[i, j, k] /= self.numped[i, j, k]
                        self.rmsped[i, j, k] /= self.numped[i, j, k]
                        self.rmsped[i, j, k] = np.sqrt(self.rmsped[i, j, k] - self.meanped[i, j, k]*self.meanped[i, j, k])
                    else:
                        print("Errot !!! Pix", i << " gain: ", j, ", capacitor ", k, " not enough events")


def remove_ped(evt, ped):
    size4drs = 4 * 1024
    for i in range(0, 8):
        for j in range(0, 2):
            fc = int(evt.firstcap[i, j])
            for k in range(0, 40): # should be roisize
                evt.samples[i, j, k] -= ped.meanped[i, j, (k+fc)%size4drs]

try:
    f = open("Randome7kHz20kev_run1.dat", "rb")
    ev = Event(40)
    ev.read(f)
    ped = PedestalSimple()
    ped.fill_ped_event(ev)
    print(ev.samples[0][1][:])
    remove_ped(ev, ped)
    print(ev.samples[0][1][:])
finally:
    f.close()