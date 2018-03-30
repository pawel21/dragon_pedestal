import numpy as np
import warnings


class PedestalSimple:
    "Calculate pedestal as a function of an absolute position in DRS4"
    n_channels = 8
    RoI = 40
    size4drs = 4*1024

    def __init__(self):
        self.channel_pedestal_value_high_gain = np.zeros((self.n_channels, self.size4drs))
        self.number_of_event_high_gain = np.zeros((self.n_channels, self.size4drs))
        self.mean_value_channel_high_gain = np.zeros((self.n_channels, self.size4drs))

    def fill_pedestal_event(self, event):
        for j in range(0, self.n_channels):
            for k in range(2, self.RoI - 2):
                # for high gain
                position = int((k + event.first_capacitor_high_gain[j]) % self.size4drs)
                self.channel_pedestal_value_high_gain[j, position] += event.samples_high_gain[j, k]
                self.number_of_event_high_gain[j, position] += 1

    def finalize_pedestal(self):
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try:
                self.mean_value_channel_high_gain = self.channel_pedestal_value_high_gain[:, :]/self.number_of_event_high_gain[:, :]
            except Warning as e:
                print("Not enough events. Error: ", e)


def remove_pedestal(event, pedestal):
    n_channels = 8
    size4drs = 4 * 1024
    for j in range(0, n_channels):
        for k in range(0, event.RoI):
            position = int((k + event.first_capacitor_high_gain[j]) % size4drs)
            event.samples_high_gain[j, k] -= pedestal.mean_value_channel_high_gain[j, position]