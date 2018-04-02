import matplotlib.pyplot as plt
import numpy as np
from protozfits.simple import File
from pedestal import PedestalSimple, remove_pedestal
from event import Event

plt.rcParams.update({'font.size': 18})


def plot_hist(event_data, gain):
    position = [(0,0), (0, 1), (0,2), (0,3),
                (1,0), (1, 1), (1,2), (1,3)]
    fig, ax = plt.subplots(2,4)
    for i in range(0, 7):
        ax[position[i]].hist(event_data[:, i, 2:38].ravel(), bins=50, facecolor='green', alpha=0.75, histtype='stepfilled')
        ax[position[i]].set_xlabel("signal [counts]")
        ax[position[i]].set_ylabel("number of events")
        sigma = np.std(event_data[:, i, 2:38].ravel())
        mu = np.mean(event_data[:, i, 2:38].ravel())
        ax[position[i]].set_title("{} Channel {}".format(gain, i + 1))
        textstr = '$\mu=%.2f$\n $\sigma=%.2f$' % (mu, sigma)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax[position[i]].text(0.05, 0.95, textstr, transform=ax[position[i]].transAxes, fontsize=14, verticalalignment='top', bbox=props)
    ax[position[7]].axis('off')
    plt.subplots_adjust(top=0.925, bottom=0.08, left=0.080, right=0.95, hspace=0.305, wspace=0.440)


def _get_pedestal(file1, N1):
    f1 = File(file1)
    ped = PedestalSimple()
    for i in range(0, N1):
        ev = next(f1.Events)
        Ev = Event(ev)
        Ev.read()
        ped.fill_pedestal_event(Ev)

    ped.finalize_pedestal()

    return ped


def subtract_pedestal(file1, file2, N1=10000, N2=3000):
    ped = _get_pedestal(file1, N1)
    f2 = File(file2)
    event_hi_gain_before_remove_pedestal = np.zeros((N2, 8, 40))
    event_hi_gain_after_remove_pedestal = np.zeros((N2, 8, 40))
    event_low_gain_before_remove_pedestal = np.zeros((N2, 8, 40))
    event_low_gain_after_remove_pedestal = np.zeros((N2, 8, 40))
    for i in range(0, N2):
        ev = next(f2.Events)
        Ev = Event(ev)
        Ev.read()
        event_hi_gain_before_remove_pedestal[i, :, :] = Ev.samples_high_gain
        event_low_gain_before_remove_pedestal[i, :, :] = Ev.samples_low_gain
        remove_pedestal(Ev, ped)
        event_hi_gain_after_remove_pedestal[i, :, :] = Ev.samples_high_gain
        event_low_gain_after_remove_pedestal[i, :, :] = Ev.samples_low_gain

    plot_hist(event_hi_gain_before_remove_pedestal, "Hi Gain")
    plot_hist(event_hi_gain_after_remove_pedestal, "Hi Gain")
    plot_hist(event_low_gain_before_remove_pedestal, "Low Gain")
    plot_hist(event_low_gain_after_remove_pedestal, "Low Gain")

subtract_pedestal("../data/Run021.1.fits.fz", "../data/Run026.1.fits.fz")
plt.show()

