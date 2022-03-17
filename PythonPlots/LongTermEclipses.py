# Preamble
import matplotlib.pyplot as plt
from math import pi
from PlotFunctions import open_csv, dark_figure, finish_figure, get_plot_path

target_folder = "LongTermEclipses"
data_path = get_plot_path(target_folder)

stela_datevector = open_csv(data_path, "stela_datevector.csv")
stela_eclipses = open_csv(data_path, "stela_eclipses.csv")
stela_mean_kep = open_csv(data_path, "stela_mean_kep.csv")

# Figures
F1, _ = dark_figure()
plt.plot(stela_datevector - stela_datevector[0], stela_eclipses / 60)
plt.xlim([0, max(stela_datevector - stela_datevector[0])])
plt.ylim([0, max(stela_eclipses) * 1.1 / 60])
plt.xlabel("Days since beginning of mission")
plt.ylabel("Eclipse duration per orbit [min]")
plt.title("Evolution of the eclipse duration over the entire mission")
finish_figure(F1, target_folder + '/' + 'long_term_eclipses.png', show=True)

F2, axes = dark_figure(subplots=(2, 3), figsize=(10, 5))
axes[0].plot(stela_datevector - stela_datevector[0], stela_mean_kep[:, 0])
axes[0].set(ylabel="SMA")
axes[1].plot(stela_datevector - stela_datevector[0], stela_mean_kep[:, 2] * 180 / pi)
axes[1].set(ylabel="Inclination [deg]")
axes[2].plot(stela_datevector - stela_datevector[0], stela_mean_kep[:, 1])
axes[2].set(ylabel="Eccentricity [-]")
axes[3].plot(stela_datevector - stela_datevector[0], stela_mean_kep[:, 3] * 180 / pi)
axes[3].set(ylabel="Argument of perigee [deg]")
axes[4].plot(stela_datevector - stela_datevector[0], stela_mean_kep[:, 4] * 180 / pi)
axes[4].set(ylabel="RAAN [deg]")
axes[5].plot(stela_datevector - stela_datevector[0], stela_mean_kep[:, 5] * 180 / pi)
axes[5].set(ylabel="Mean anomaly [deg]")
plt.suptitle("Evolution of orbital parameters over the entire mission", color='white')
finish_figure(F2, target_folder + '/' + 'long_term_orbit.png', show=True)
