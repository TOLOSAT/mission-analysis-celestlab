# Preamble
import matplotlib.pyplot as plt
from math import pi
import numpy as np
from PlotFunctions import open_csv, dark_figure, finish_figure, get_plot_path
import cartopy.crs as ccrs

target_folder = "PrettyVisuals"
data_path = get_plot_path(target_folder)

cjd = np.array(open_csv(data_path, "lydane_cjd.csv")).transpose()
cjd0 = np.array(open_csv(data_path, "cjd0.csv")).transpose()
ecc = np.array(open_csv(data_path, "lydane_ecc.csv")).transpose()
inc = np.array(open_csv(data_path, "lydane_inc.csv")).transpose()
pom = np.array(open_csv(data_path, "lydane_pom.csv")).transpose()
RAAN = np.array(open_csv(data_path, "lydane_RAAN.csv")).transpose()
sma = np.array(open_csv(data_path, "lydane_sma.csv")).transpose()
meanAnomaly = np.array(open_csv(data_path, "lydane_M.csv")).transpose()
pos_lla = np.array(open_csv(data_path, "pos_lla.csv")).transpose()
pos_lla[:, ::1] = pos_lla[:, ::1] * 180 / pi

earthRadius = 6.3781e6

# Figures
F1, axes = dark_figure(figsize=(7, 4.5))
axes[0].axis('off')
axes = plt.axes(projection=ccrs.PlateCarree())
axes.stock_img()
axes.set_global()
axes.plot(pos_lla[:, 0], pos_lla[:, 1], color='blue', linewidth=1, transform=ccrs.Geodetic())
# gl = axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                   xlocs=np.arange(-150,360,30), ylocs=np.arange(-30, 40,10),
#                   x_inline=False, y_inline=False, linewidth=0.33, color='k',alpha=0.5)
# gl.right_labels = False
# gl.top_labels = False
plt.suptitle("Ground track across " + str(int(cjd[-1][0] - cjd0[0])) + " day(s)", color='white')
finish_figure(F1, target_folder + '/' + target_folder + '_GroundTrack.png', show=True)

tmp_xTime = (cjd - cjd0)
F2, axes = dark_figure(subplots=(2, 3), figsize=(10, 6.5))
handles = []
axes[0].plot(tmp_xTime, (sma - earthRadius) / 1e3)
axes[1].plot(tmp_xTime, inc * 180 / pi)
axes[2].plot(tmp_xTime, ecc)
axes[3].plot(tmp_xTime, pom * 180 / pi)
axes[4].plot(tmp_xTime, RAAN * 180 / pi)
axes[5].plot(tmp_xTime, meanAnomaly * 180 / pi, label='dummy label')
for jj in range(6):
    axes[jj].set(xlabel="Time [days]", xlim=[0, axes[jj].get_xlim()[1]])
axes[0].set(ylabel="Altitude [km]")
axes[1].set(ylabel="Inclination [deg]")
axes[2].set(ylabel="Eccentricity [-]")
axes[3].set(ylabel="Argument of perigee [deg]")
axes[4].set(ylabel="RAAN [deg]")
axes[5].set(ylabel="Mean Anomaly [deg]")
plt.suptitle("Evolution of orbital parameters", color='white')
finish_figure(F2, target_folder + '/' + target_folder + '_OrbitalParameters.png', show=True)
print("done")
