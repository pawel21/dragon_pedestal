import matplotlib.pyplot as plt
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
                        print("error pix", self.numped[i, j, k])
                        #print("Errot !!! Pix", i << " gain: ", j, ", capacitor ", k, " not enough events")


def remove_ped(evt, ped):
    size4drs = 4 * 1024
    for i in range(0, 8):
        for j in range(0, 2):
            fc = int(evt.firstcap[i, j])
            for k in range(0, 40): # should be roisize
                evt.samples[i, j, k] -= ped.meanped[i, j, (k+fc)%size4drs]

try:
    f1 = open("Randome7kHz20kev_run1.dat", "rb")
    ev1 = Event(40)
    ped = PedestalSimple()
    for i in range(1, 2000):
        f1.seek((i*1344))
        ev1.read(f1)
        ped.fill_ped_event(ev1)
    ped.finalize_ped()
    f2 = open("Randome7kHz20kev_run2.dat", "rb")
    ev2 = Event(40)
    for i in range(1, 20):
        f2.seek((i * 1344))
        ev2.read(f2)
    ev2_before_remove_ped = ev2.samples
    ev_before = ev2.samples[0][0][:]
    fig, ax = plt.subplots(2)
    ax[0].plot(ev_before)
    print(ev_before)
    remove_ped(ev2, ped)
    ev_after = ev2.samples[0][0][:]
    print(ev_after)
    ev2_after_remove_ped = ev2.samples
    ax[1].plot(ev_after)
    plt.show()

    #statistics
    nn = 7
    window = 6
    hped = np.zeros((nn, 2))
    for i in range(0, nn):
        for j in range(0, 2):
            for k in range(0, 31):
                hped[i][j] = ev2.SumSlices(i, j, k, window)/(np.sqrt(1.*window))
finally:
    f1.close()
    f2.close()