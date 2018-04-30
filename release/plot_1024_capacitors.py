import matplotlib.pyplot as plt
import numpy as np
import os
from protozfits import SimpleFile

from event import Event
from pedestal import PedestalSimple

plt.rcParams.update({'font.size': 35})

f = SimpleFile(os.path.join("..", "data", "Run021.1.fits.fz"))
N = 10000

ped = PedestalSimple()

for i in range(1, N):
    ev = next(f.Events)
    Ev = Event(ev)
    Ev.read()
    ped.fill_pedestal_event(Ev)

ped.finalize_pedestal()

x = np.linspace(1, 1024, 1024)

plt.errorbar(x, ped.mean_value_channel_high_gain[0, :1024], yerr=ped.rms_high_gain[0, :1024], fmt='o')
plt.xlabel("DRS 4 kondensator")
plt.ylabel("Offset komórkki wraz z RMS")
plt.grid(True)
plt.show()