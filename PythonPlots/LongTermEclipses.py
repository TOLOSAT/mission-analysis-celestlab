# Preamble
from datetime import datetime, timezone
import matplotlib.pyplot as plt
import numpy as np
import os
import csv
from PIL import Image

Badge_TOLOSAT = Image.open('assets/TOLOSAT.png')
target_folder = "LongTermEclipses"
file_path = os.path.dirname(os.path.realpath(__file__)).replace("PythonPlots", target_folder)
os.makedirs("../" + target_folder, exist_ok=True)


def open_csv(path, target_file):
    tmp = []
    with open(path + "\\" + target_file, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if len(row) == 1:
                tmp.append(row[0])
            else:
                tmp.append(row)
    return tmp


stela_datevector = open_csv(file_path, "stela_datevector.csv")
stela_eclipses = open_csv(file_path, "stela_eclipses.csv")
stela_mean_kep = open_csv(file_path, "stela_mean_kep.csv")


# Functions for figures
def dark_figure():
    fig = plt.figure(facecolor='#0D1117', figsize=(7, 5.2))
    f = plt.gca()
    for i in f.spines:
        f.spines[i].set_color('white')
    f.tick_params(axis='x', colors='white', which='both')
    f.tick_params(axis='y', colors='white', which='both')
    f.yaxis.label.set_color('white')
    f.xaxis.label.set_color('white')
    f.title.set_color('white')
    f.set_facecolor('#0D1117')
    return f, fig


def finish_figure(fig, path, show):
    plt.tight_layout()
    fig.subplots_adjust(bottom=0.20)
    fig_axes1 = fig.add_axes([0.678, 0.02, 0.3, 0.3], anchor='SE', zorder=1)
    fig_axes1.imshow(Badge_TOLOSAT)
    fig_axes1.axis('off')
    plt.savefig('plots/' + path + '_transparent.png', transparent=True, dpi=500)
    plt.savefig('plots/' + path + '_background.png', transparent=False, dpi=500)
    if show:
        plt.show()
    plt.close()


def flip_legend(reverse):
    handles_, labels_ = plt.gca().get_legend_handles_labels()
    handles_ = [k for j in [handles_[i::4] for i in range(4)] for k in j]
    labels_ = [k for j in [labels_[i::4] for i in range(4)] for k in j]
    if reverse:
        return handles_[::-1], labels_[::-1]
    else:
        return handles_, labels_


def flatten(list_of_lists):
    flattened_list = []
    for i in list_of_lists:
        if isinstance(i, list):
            flattened_list += i
        else:
            flattened_list.append(i)
    return flattened_list


# Figures
_, F1 = dark_figure()
plt.plot(stela_datevector, stela_eclipses)
plt.xlabel(datetime.now(timezone.utc).strftime("Plot generated on %Y/%m/%d at %H:%M:%S UTC."), color='dimgray',
           labelpad=10)
finish_figure(F1, target_folder + 'test1.png', show=True)
