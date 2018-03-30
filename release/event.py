import numpy as np
from ctapipe.core import Container, Field

from tools import translate_fits


class ContainerEvent(Container):
    n_channels = 8
    RoI = 40
    hi_gain_samples = Field(np.zeros((n_channels, RoI)), "Hi Gain Samples")


class Event:
    n_channels = 8
    RoI = 40

    def __init__(self, event):
        self.event = event
        self.samples_high_gain = np.zeros((self.n_channels, self.RoI))
        self.samples_low_gain = np.zeros((self.n_channels, self.RoI))
        self.first_capacitor_high_gain = np.zeros(8)
        self.first_capacitor_low_gain = np.zeros(8)

    def read(self):
        hi_gain, low_gain = translate_fits(self.event)
        for j in range(0, self.n_channels):
            self.samples_high_gain[j, 0:40] = hi_gain[j * 40:(j + 1) * 40]
            self.samples_low_gain[j, 0:40] = low_gain[j * 40:(j + 1) * 40]
        fc = self.event.firstCapacitorIds.data[0:8]
        for i in [0, 2, 4, 6]:
            self.first_capacitor_high_gain[i] = fc[i]
            self.first_capacitor_high_gain[i + 1] = fc[i]
        for j in [1, 3, 5, 7]:
            self.first_capacitor_low_gain[j - 1] = fc[j]
            self.first_capacitor_low_gain[j] = fc[j]





