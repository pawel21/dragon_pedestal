import matplotlib.pyplot as plt
import numpy as np
from protozfits import SimpleFile
from pedestal import PedestalSimple, remove_pedestal
from event import Event
from tools import plot_hist

plt.rcParams.update({'font.size': 32})

N1 = 15000
f1 = SimpleFile("../data/Run021.1.fits.fz")
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
    fc_hg_old = Ev.first_capacitor_high_gain
    ev = next(f1.Events)
    Ev = Event(ev)
    Ev.read()
    remove_pedestal(Ev, ped)
    event_hi_gain_after_remove_pedestal[i, :, :] = Ev.samples_high_gain

    for j in range(0, 8):
        # looking for spike A
        for k in range(0, 4):
            abspos = int(roisize-2 + fc_hg_old[j]+ k*1024)
            pos = int((abspos - Ev.first_capacitor_high_gain[j] + size4drs)%size4drs)
            if ((pos > 0) and (pos <= roisize-1)) or (pos == size4drs -1):
                print("find something !!!!!!!!!!!")
                print(pos)
                t = np.arange(0, 40, 1)
                plt.step(t, event_hi_gain_after_remove_pedestal[i, j, :], 'g-')
                plt.ylabel("signal [counts]")
                plt.xlabel("time sample [ns]")
                plt.title("Spike type A")
                plt.grid(True)
                plt.show()

            abspos = int(1024 - roisize-2 - fc_hg_old[j]+ k*1024 +size4drs)
            pos = int((abspos - Ev.first_capacitor_high_gain[j] + size4drs)%size4drs)
            if ((pos > 0) and (pos <= roisize-1)) or (pos == size4drs -1):
                print("find something")
                print(pos)

                if pos>2 and pos<38:
                    t = np.arange(2, 38, 1)
                    samples = event_hi_gain_after_remove_pedestal[i, j, 2:38]
                    pos = pos - 2
                    fig, ax1 = plt.subplots(2, 1)
                    ax1[0].step(t, samples, 'g-', lw=3)
                    ax1[0].set_ylabel("sygnał")
                    ax1[0].set_xlabel("czas [ns]")
                    ax1[0].grid(True)
                    ax1[0].set_ylim([-20, 80])

                    samples[pos] = samples[pos - 1] + 0.33*(samples[pos + 2] - samples[pos-1])
                    samples[pos+1] = samples[pos -1] + 0.66*(samples[pos + 2] - samples[pos - 1])

                    ax1[1].step(t, samples, 'g-', lw=3)
                    ax1[1].set_ylabel("sygnał")
                    ax1[1].set_xlabel("czas [ns]")
                    ax1[1].grid(True)
                    ax1[1].set_ylim([-20, 80])

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
            t = np.arange(2, 38, 1)

            samples = event_hi_gain_after_remove_pedestal[i, j, 2:38]


            fig, ax1 = plt.subplots(2, 1)
            ax1[0].step(t, samples, 'g-', lw=3)
            ax1[0].set_ylabel("sygnał")
            ax1[0].set_xlabel("czas [ns]")
            ax1[0].grid(True)
            ax1[0].set_ylim([-20, 40])


            samples[spike_b_pos] = 0.5*(samples[spike_b_pos-1] + samples[spike_b_pos+1])

            ax1[1].step(t, samples, 'g-', lw=3)
            ax1[1].set_ylabel("sygnał")
            ax1[1].set_xlabel("czas [ns]")
            ax1[1].grid(True)
            ax1[1].set_ylim([-20, 40])


            plt.show()

    plt.close("all")
plot_hist(event_hi_gain_before_remove_pedestal, "Before HG")
plot_hist(event_hi_gain_after_remove_pedestal, "After HG")
plt.show()