# Preamble
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from math import pi
import numpy as np
from PlotFunctions import open_csv, dark_figure, finish_figure, get_plot_path, flatten, flip_legend

target_folder = "LaunchYearMissionDuration"
data_path = get_plot_path(target_folder)

cjd = np.array(open_csv(data_path, "cjd_stela.csv")).transpose()
cjd0 = np.array(open_csv(data_path, "cjd0.csv")).transpose()
ecc = np.array(open_csv(data_path, "ecc.csv")).transpose()
inc = np.array(open_csv(data_path, "inc.csv")).transpose()
mltan = np.array(open_csv(data_path, "mltan.csv")).transpose()
pom = np.array(open_csv(data_path, "pom.csv")).transpose()
RAAN = np.array(open_csv(data_path, "RAAN.csv")).transpose()
sma = np.array(open_csv(data_path, "sma.csv")).transpose()
years = np.array(open_csv(data_path, "years.csv")).transpose()
earthRadius = 6.3781e6

# colors = np.array([(0, 135, 108), (61, 154, 112), (100, 173, 115), (137, 191, 119), (175, 209, 124), (214, 225, 132),
#                    (255, 241, 143), (253, 213, 118), (251, 184, 98), (245, 155, 86), (238, 125, 79), (227, 94, 78),
#                    (212, 61, 81)]) / 255.0
colors = pl.cm.jet(np.linspace(0, 1, len(years)))

# Figures
F1, axes = dark_figure(subplots=(2, 3), figsize=(10, 6))
handles = []
for ii in range(len(years)):
    tmp_xTime = cjd[:, ii] - cjd0[ii]
    color = colors[ii]
    axes[0].plot(tmp_xTime, (sma[:, ii] - earthRadius) / 1e3, color=color)
    axes[1].plot(tmp_xTime, inc[:, ii] * 180 / pi, color=color)
    axes[2].plot(tmp_xTime, ecc[:, ii], color=color)
    axes[3].plot(tmp_xTime, pom[:, ii] * 180 / pi, color=color)
    axes[4].plot(tmp_xTime, RAAN[:, ii] * 180 / pi, color=color)
    axes[5].plot(tmp_xTime, mltan[:, ii], color=color, label='dummy label')
axes[0].set(ylabel="Altitude [km]")
axes[1].set(ylabel="Inclination [deg]")
axes[2].set(ylabel="Eccentricity [-]")
axes[3].set(ylabel="Argument of perigee [deg]")
axes[4].set(ylabel="RAAN [deg]")
axes[5].set(ylabel="MLTAN [hours]")
plt.suptitle("Evolution of orbital parameters over the entire mission for various launch years", color='white')
handles, _ = axes[5].get_legend_handles_labels()
handles, labels = flip_legend(7, False, handles, [str(int(x)) for x in flatten(years.tolist())])
F1.legend(handles, labels, loc=(0.015, 0.055), ncol=7, frameon=False,
          labelcolor='white')
finish_figure(F1, target_folder + '/' + 'launchYearMissionDuration.png', show=True)
print("done")
