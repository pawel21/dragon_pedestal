import matplotlib.pyplot as plt
import numpy as np
from protozfits import SimpleFile
from pedestal import PedestalSimple, remove_pedestal
from event import Event
from tools import plot_hist

plt.rcParams.update({'font.size': 32})

N1 = 15000
f1 = SimpleFile("../data/Run027.1.fits.fz")
ped = PedestalSimple()

for i in range(0, N1):
    ev = next(f1.Events)
    Ev = Event(ev)
    Ev.read()
    ped.fill_pedestal_event(Ev)

ped.finalize_pedestal()


N2 = 5000
roisize = 40
size4drs = 4*1024

event_hi_gain_before_remove_pedestal = np.zeros((N2, 8, 40))
event_hi_gain_after_remove_pedestal = np.zeros((N2, 8, 40))
fc_hg_old = np.zeros(8)

ev = next(f1.Events)
Ev = Event(ev)
Ev.read()

for i in range(0, N2):
    event_hi_gain_before_remove_pedestal[i, :, :] = Ev.samples_high_gain
    remove_pedestal(Ev, ped)
    event_hi_gain_after_remove_pedestal[i, :, :] = Ev.samples_high_gain
    fc_hg_old = Ev.first_capacitor_high_gain
    ev = next(f1.Events)
    Ev = Event(ev)
    Ev.read()
    for j in range(0, 8):

        # looking for spike A
        for k in range(0, 4):
            abspos = int(roisize-2 + fc_hg_old[j]+k*1024)
            pos = (abspos - Ev.first_capacitor_high_gain[j] + size4drs)%size4drs
            if ((pos > 0) and (pos <= roisize-1)) or (pos == size4drs -1):
                print("find something !!!!!!!!!!!")
                t = np.arange(2, 38, 1)
                plt.step(t, event_hi_gain_after_remove_pedestal[i, j, 2:38], 'g-')
                plt.ylabel("signal [counts]")
                plt.xlabel("time sample [ns]")
                plt.title("Spike type A")
                plt.grid(True)
                plt.show()

            abspos = int(1024 - roisize-2 -fc_hg_old[j]+k*1024+size4drs)
            pos = int((abspos - Ev.first_capacitor_high_gain[j] + size4drs)%size4drs)
            if ((pos > 0) and (pos <= roisize-1)) or (pos == size4drs -1):
                print("find something !!!!!!!!!!!")
                plt.step(event_hi_gain_after_remove_pedestal[i, j, 2:38], 'g-')
                plt.ylabel("signal [counts]")
                plt.xlabel("time sample [ns]")
                plt.title("Spike type A")
                plt.grid(True)
                plt.show()

        spike_b_pos = int((fc_hg_old[j] - 1 - Ev.first_capacitor_high_gain[j] + 2*size4drs)%size4drs)
        print(spike_b_pos)
        if spike_b_pos == 0:
            plt.step(event_hi_gain_after_remove_pedestal[i, j, 2:38], 'g-')
            plt.ylabel("signal [counts]")
            plt.xlabel("time sample [ns]")
            plt.grid(True)
            plt.ylim([-20, 110])
            plt.title("Spike type B")
            plt.show()
        if spike_b_pos == roisize-1:
            plt.step(event_hi_gain_after_remove_pedestal[i, j, 2:38], 'g-')
            plt.ylabel("signal [counts]")
            plt.xlabel("time sample [ns]")
            plt.grid(True)
            plt.ylim([-20, 110])
            plt.title("Spike type B")
            plt.show()
        if spike_b_pos < roisize-1:
            plt.step(event_hi_gain_after_remove_pedestal[i, j, 2:38], 'g-')
            plt.ylabel("signal [counts]")
            plt.xlabel("time sample [ns]")
            plt.grid(True)
            plt.ylim([-20, 110])
            plt.title("Spike type B")
            plt.show()

plot_hist(event_hi_gain_before_remove_pedestal, "Before HG")
plot_hist(event_hi_gain_after_remove_pedestal, "After HG")
plt.show()