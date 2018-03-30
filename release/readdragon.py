import matplotlib.pyplot as plt
import numpy as np
from protozfits.simple import File
from pedestal import PedestalSimple, remove_pedestal
from event import Event

plt.rcParams.update({'font.size': 20})

f = File("../data/Run021.1.fits.fz")
N = 10000
ped = PedestalSimple()

for i in range(0, N):
    ev = next(f.Events)
    Ev = Event(ev)
    Ev.read()
    ped.fill_pedestal_event(Ev)

ped.finalize_pedestal()


f2 = File("../data/Run027.2.fits.fz")
N2 = 2000
event_hi_gain_before_remove_pedestal = np.zeros((N2, 8, 40))
event_hi_gain_after_remove_pedestal = np.zeros((N2, 8, 40))
for i in range(0, N2):
    ev = next(f2.Events)
    Ev = Event(ev)
    Ev.read()
    event_hi_gain_before_remove_pedestal[i, :, :] = Ev.samples_high_gain[:, :]
    remove_pedestal(Ev, ped)
    event_hi_gain_after_remove_pedestal[i, :, :] = Ev.samples_high_gain




position = [(0,0), (0, 1), (0,2), (0,3),
            (1,0), (1, 1), (1,2), (1,3)]

fig, ax = plt.subplots(2,4)
for i in range(0, 8):
    ax[position[i]].hist(event_hi_gain_before_remove_pedestal[:, i, 2:38].ravel(), bins=50, facecolor='green', alpha=0.75)
    sigma = np.std(event_hi_gain_before_remove_pedestal[:, i, 2:38].ravel().flatten())
    mu = np.mean(event_hi_gain_before_remove_pedestal[:, i, 2:38].ravel().flatten())
    ax[position[i]].set_title("Hi Gain Channel {}".format((i + 1)))
    textstr = '$\mu=%.2f$\n$\sigma=%.2f$' % (mu, sigma)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax[position[i]].text(0.05, 0.95, textstr, transform=ax[position[i]].transAxes, fontsize=14, verticalalignment='top', bbox=props)



plt.figure()
for i in range(0, 8):
    plt.subplot(2, 4, i+1)
    plt.hist(event_hi_gain_before_remove_pedestal[:, i, 2:38].ravel(), bins=50, facecolor='green', alpha=0.75)
    plt.xlabel("signal [counts]")
    plt.ylabel("number of events")
    sigma = np.std(event_hi_gain_before_remove_pedestal[:, i, 2:38].ravel().flatten())
    mu = np.mean(event_hi_gain_before_remove_pedestal[:, i, 2:38].ravel().flatten())
    textstr = '$\mu=%.2f$\n$\sigma=%.2f$' % (mu, sigma)
    plt.title("Hi Gain Channel {}".format((i + 1)))

plt.legend()

plt.figure()
for i in range(0, 8):
    plt.subplot(2, 4, i+1)
    plt.hist(event_hi_gain_after_remove_pedestal[:, i, 2:38].ravel(), bins=50, facecolor='green', alpha=0.75)
    plt.xlabel("signal [counts]")
    plt.ylabel("number of events")
    sigma = np.std(event_hi_gain_after_remove_pedestal[:, i, 2:38].ravel().flatten())
    mu = np.mean(event_hi_gain_after_remove_pedestal[:, i, 2:38].ravel().flatten())
    textstr = '$\mu=%.2f$\n$\sigma=%.2f$' % (mu, sigma)
    plt.title("Hi Gain Channel {}".format((i + 1)))
plt.legend()

plt.show()