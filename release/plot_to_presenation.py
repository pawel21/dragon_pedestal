import matplotlib.pyplot as plt
import numpy as np
from protozfits.simple import File
from pedestal import PedestalSimple, remove_pedestal
from event import Event

plt.rcParams.update({'font.size': 32})

N1 = 20000
f1 = File("../data/Run021.1.fits.fz")
ped = PedestalSimple()

for i in range(0, N1):
    ev = next(f1.Events)
    Ev = Event(ev)
    Ev.read()
    ped.fill_pedestal_event(Ev)

ped.finalize_pedestal()

N2 = 1000
event_hi_gain_before_remove_pedestal = np.zeros((N2, 8, 40))
event_hi_gain_after_remove_pedestal = np.zeros((N2, 8, 40))


for i in range(0, N2):
    ev = next(f1.Events)
    Ev = Event(ev)
    Ev.read()
    event_hi_gain_before_remove_pedestal[i, :, :] = Ev.samples_high_gain
    remove_pedestal(Ev, ped)
    event_hi_gain_after_remove_pedestal[i, :, :] = Ev.samples_high_gain

fig, (ax0, ax1) = plt.subplots(ncols = 2)

sigma_before = np.std(event_hi_gain_before_remove_pedestal[:, 0, 2:38].ravel())
mu_before = np.mean(event_hi_gain_before_remove_pedestal[:, 0, 2:38].ravel())
ax0.hist(event_hi_gain_before_remove_pedestal[:, 0, 2:38].ravel(), bins=50)
ax0.set_xlabel("signal [counts]")
ax0.set_ylabel("number of events")
ax0.set_title("Before")
textstr = '$\mu=%.2f$\n $\sigma=%.2f$' % (mu_before, sigma_before)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax0.text(0.05, 0.95, textstr, transform=ax0.transAxes, fontsize=32, verticalalignment='top', bbox=props)
ax0.set_xlim([195, 405])
ax0.grid()

sigma_after = np.std(event_hi_gain_after_remove_pedestal[:, 0, 2:38].ravel())
mu_after = np.mean(event_hi_gain_after_remove_pedestal[:, 0, 2:38].ravel())
ax1.hist(event_hi_gain_after_remove_pedestal[:, 0, 2:38].ravel(), bins=50)
ax1.set_xlabel("signal [counts]")
ax1.set_ylabel("number of events")
ax1.set_title("After")
textstr = '$\mu=%.2f$\n $\sigma=%.2f$' % (mu_after, sigma_after)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax1.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=32, verticalalignment='top', bbox=props)
ax1.set_xlim([-105, 105])
ax1.grid()

fig, (ax0, ax1) = plt.subplots(nrows = 2)

ax0.plot(event_hi_gain_before_remove_pedestal[3, 0, 0:40], 'b-', lw=3)
ax0.set_ylabel("signal [counts]")
ax0.set_xlabel("time sample [ns]")
ax0.axhline(color='k', lw=1)
ax0.set_ylim([-80, 350])
ax0.grid()

ax1.plot(event_hi_gain_after_remove_pedestal[3, 0, 0:40], 'b-', lw=3)
ax1.set_ylabel("signal [counts]")
ax1.set_xlabel("time sample [ns]")
ax1.set_ylim([-80, 350])
ax1.axhline(color='k', lw=1)
plt.grid(True)
plt.show()
