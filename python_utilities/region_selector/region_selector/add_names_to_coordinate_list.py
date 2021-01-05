"""
Ad-hoc adding of arbitrary labels to regions produced using region_selector.py.
Shows RGB images with overlays over regions and associated labels.
"""

import sys
import pickle
from region_selector import get_show_image
import matplotlib.pyplot as plt
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Show/set label on points/regions.')
parser.add_argument('file', metavar='file', type=str, help='Path to region file.')
parser.add_argument('--new-labels', type=str, default=None, help='List of new labels, on format label1,label2,label3,... . Regions are drawn in order, and labels drawn from this list. If not specified, will just use already contained labels.')
parser.add_argument('--no-overwrite', default=False, action='store_true', help="Don't write to point file again.")
parser.add_argument('--legend', default=False, action='store_true')
args = parser.parse_args()

plt.ioff()

coordlist = pickle.load(open(args.file, 'rb'))

#get labels from cli or use already collected labels
if args.new_labels is not None:
    sample_names = args.new_labels.split(',')
else:
    try:
        sample_names = list(np.concatenate([entry['labels'] for entry in coordlist]))
    except:
        sample_names = []

import matplotlib.patches as patches

for i,entry in enumerate(coordlist):
    fname = entry['filename']

    #show image
    plt.figure(1); plt.clf()
    plt.imshow(get_show_image(fname))

    #plot points and add labels
    labels = []

    is_coords = True
    length = 0
    try:
        coordinates = entry['coordinates']
        length = len(coordinates)
    except:
        is_coords = False
        pos = entry['skin_pos']
        size = entry['skin_size']
        length = len(pos)

    for j in range(length):
        try:
            label = sample_names.pop(0)
        except:
            label=None

        if is_coords:
            coord = coordinates[j]
            plt.plot(coord[0], coord[1], 'o', label=label)
            plt.text(coord[0], coord[1], label)
        else:
            plt.gca().add_patch(patches.Rectangle(pos[j], size[j][0], size[j][1], facecolor='none', edgecolor='black'))
            plt.text(pos[j][0], pos[j][1], label)

        labels.append(label)

    if args.legend:
        plt.legend()
    plt.show()

    coordlist[i]['labels'] = labels

#save again with new labels (if any)
if not args.no_overwrite:
    pickle.dump(coordlist, open(args.file, "wb"))
