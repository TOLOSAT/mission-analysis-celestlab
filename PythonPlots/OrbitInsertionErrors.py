# Preamble
import matplotlib.pyplot as plt
from math import pi
import numpy as np
from PlotFunctions import open_csv, dark_figure, finish_figure, get_plot_path, flip_legend

target_folder = "OrbitInsertionErrors"
data_path = get_plot_path(target_folder)

cjd = np.array(open_csv(data_path, "cjd_stela.csv")).transpose()
cjd0 = np.array(open_csv(data_path, "cjd0.csv")).transpose()
ecc = np.array(open_csv(data_path, "ecc_stela.csv")).transpose()
inc = np.array(open_csv(data_path, "inc_stela.csv")).transpose()
mltan = np.array(open_csv(data_path, "mltan.csv")).transpose()
pom = np.array(open_csv(data_path, "pom_stela.csv")).transpose()
RAAN = np.array(open_csv(data_path, "RAAN_stela.csv")).transpose()
sma = np.array(open_csv(data_path, "sma_stela.csv")).transpose()
eclipses = np.array(open_csv(data_path, "eclipse_duration_umb.csv")).transpose()
labels = open_csv(data_path, "legend_name.csv", is_string=True)
labels = [ii.replace('Â', '') for ii in labels]

earthRadius = 6.3781e6

# Figures
tmp_xTime = (cjd - cjd0) / 365.25

F1, axes = dark_figure(figsize=(7, 6))
axes[0].plot(tmp_xTime, eclipses / 60, label='dummy label')
plt.xlim([0, axes[0].get_xlim()[1]])
plt.ylim([0, axes[0].get_ylim()[1]])
plt.xlabel("Time [years]")
plt.ylabel("Eclipse duration per orbit [min]")
plt.title("Evolution of the eclipse duration over the entire mission \n for maximum insertion errors")
handles, _ = axes[0].get_legend_handles_labels()
handles, labels = flip_legend(2, False, handles, labels)
F1.legend(handles, labels, loc=(0.015, 0.04), ncol=2, frameon=False,
          labelcolor='white')
finish_figure(F1, target_folder + '/' + target_folder + '_Eclipses.png', show=True)

F2, axes = dark_figure(subplots=(2, 3), figsize=(10, 6.5))
handles = []
for ii in range(len(labels)):
    axes[0].plot(tmp_xTime, (sma[:, ii] - earthRadius) / 1e3)
    axes[1].plot(tmp_xTime, inc[:, ii] * 180 / pi)
    axes[2].plot(tmp_xTime, ecc[:, ii])
    axes[3].plot(tmp_xTime, pom[:, ii] * 180 / pi)
    axes[4].plot(tmp_xTime, RAAN[:, ii] * 180 / pi)
    axes[5].plot(tmp_xTime, mltan[:, ii], label='dummy label')
for jj in range(6):
    axes[jj].set(xlabel="Time [years]", xlim=[0, axes[jj].get_xlim()[1]])
axes[0].set(ylabel="Altitude [km]")
axes[1].set(ylabel="Inclination [deg]")
axes[2].set(ylabel="Eccentricity [-]")
axes[3].set(ylabel="Argument of perigee [deg]")
axes[4].set(ylabel="RAAN [deg]")
axes[5].set(ylabel="MLTAN [hours]")
plt.suptitle("Evolution of orbital parameters during the entire mission for maximum insertion errors", color='white')
handles, _ = axes[5].get_legend_handles_labels()
handles, labels = flip_legend(2, False, handles, labels)
F2.legend(handles, labels, loc=(0.015, 0.055), ncol=2, frameon=False,
          labelcolor='white')
finish_figure(F2, target_folder + '/' + target_folder + '_OrbitalParameters.png', show=True)
print("done")
