import matplotlib.pyplot as plt
import numpy as np
import os
from protozfits import SimpleFile

from event import Event


def plot_hist_time_between_event(file_name, N):
    f = SimpleFile(os.path.join("..", "data", file_name))
    time = np.zeros(N)
    diff_time = np.zeros(N - 1)
    for i in range(1, N-1):
        ev = next(f.Events)
        Ev = Event(ev)
        Ev.read()
        time[i] = Ev.local133MHzClockCounter/(133*1e6)
        diff_time[i] = time[i] - time[i - 1]

    plt.hist(diff_time[2:], bins=15)
    plt.title(file_name)
    plt.show()
    return diff_time


list_of_events_file = os.listdir("../data")
for f in list_of_events_file:
    try:
        plot_hist_time_between_event(f, 5000)
    except Exception as err:
        print(err)

