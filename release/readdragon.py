import matplotlib.pyplot as plt
import numpy as np
from protozfits.simple import File
from pedestal import PedestalSimple, remove_pedestal
from event import Event


f = File("../data/Run021.1.fits.fz")
N = 40000
ped = PedestalSimple()

for i in range(0, N):
    ev = next(f.Events)
    Ev = Event(ev)
    Ev.read()
    ped.fill_pedestal_event(Ev)

ped.finalize_pedestal()


f2 = File("../data/Run027.2.fits.fz")
N2 = 5000
event_hi_gain_before_remove_pedestal = np.zeros((N2, 8, 40))
event_hi_gain_after_remove_pedestal = np.zeros((N2, 8, 40))
for i in range(0, N2):
    ev = next(f2.Events)
    Ev = Event(ev)
    Ev.read()
    event_hi_gain_before_remove_pedestal[i, :, :] = Ev.samples_high_gain[:, :]
    remove_pedestal(Ev, ped)
    event_hi_gain_after_remove_pedestal[i, :, :] = Ev.samples_high_gain


plt.figure()
for i in range(0, 8):
    plt.subplot(2, 4, i+1)
    plt.hist(event_hi_gain_before_remove_pedestal[:, i, 2:38].ravel(), bins=50, facecolor='green', alpha=0.75)
    std = np.std(event_hi_gain_before_remove_pedestal[:, i, 2:38].ravel().flatten())
    plt.title("Hi Gain Channel {}, $\sigma$ = {:.2f}.".format((i + 1), std))
plt.legend()

plt.figure()
for i in range(0, 8):
    plt.subplot(2, 4, i+1)
    plt.hist(event_hi_gain_after_remove_pedestal[:, i, 2:38].ravel(), bins=50, facecolor='green', alpha=0.75)
    std = np.std(event_hi_gain_after_remove_pedestal[:, i, 2:38].ravel())
    plt.title("Hi Gain Channel {}, $\sigma$ = {:.2f}.".format((i + 1), std))
plt.legend()

plt.show()