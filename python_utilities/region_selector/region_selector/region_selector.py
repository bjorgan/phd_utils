"""
Toolkit for selecting regions in multiple images and producing a list over
coordinates.

Originally developed for selecting skin/background regions for training a skin
classifier (therefore has labels "skin"/"background"), but later ad-hoc
extended using add_names_to_coordinate_list.py to enable arbitrary labels for
obtaining data for other classification tasks.

2019-02-14, NOTE: When loading back file list to produce segmentation
(random_forest_from_regions.py), it chooses only the first region within each
region array. Selecting multiple boxes therefore technically has no effect.
"""

import argparse
import matplotlib.pyplot as plt
from hyperspectral_utils import hyperread

import matplotlib.patches as patches
import numpy as np

from matplotlib.widgets import RectangleSelector
import imageio

import pickle

class region_selection_handler(object):
    """
    For selecting rectangles in image and collecting them.
    """

    def __init__(self):
        self.background_regions = []
        self.skin_regions = []
        self.selecting_skin = True
        self.rectangle_selector = RectangleSelector(plt.gca(), self.onselect, drawtype='box')
        plt.connect('key_press_event', self.toggle_selector)
        self.drawn_rects = []
        plt.title('Press n for skin (default), b for background, m to clear')

    def onselect(self, eclick, erelease):
        """
        Select rectangles, paint them and put in correct position array dependent on whether
        skin or background is to be selected.
        """
        start_pos = (eclick.xdata, eclick.ydata)
        end_pos = (erelease.xdata, erelease.ydata)
        rect_size = np.array(end_pos) - np.array(start_pos)

        if self.selecting_skin:
            self.skin_regions.append((start_pos, rect_size))
            color = 'yellow'
        else:
            self.background_regions.append((start_pos, rect_size))
            color = 'red'

        self.drawn_rects.append(plt.gca().add_patch(patches.Rectangle(start_pos, rect_size[0], rect_size[1], facecolor=color, alpha=0.7)))

    def clear(self):
        """
        Clear figure of rectangles, remove current rectangles.
        """
        self.background_regions.clear()
        self.skin_regions.clear()
        for rect in self.drawn_rects:
            rect.remove()
        self.drawn_rects.clear()
        plt.gcf().canvas.draw()

    def toggle_selector(self, event):
        """
        Handle key input, toggle whether skin or background is to be selected.
        """
        if event.key is 'b':
            self.selecting_skin = False
        if event.key is 'n':
            self.selecting_skin = True
        if event.key is 'm':
            self.clear()
    def to_dict(self):
        seg = {}
        seg['skin_pos'] = [region[0] for region in self.skin_regions]
        seg['skin_size'] = [region[1] for region in self.skin_regions]
        seg['background_pos'] = [region[0] for region in self.background_regions]
        seg['background_size'] = [region[1] for region in self.background_regions]
        return seg

class point_selection_handler(object):

    def __init__(self):
        self.coordinates = []
        self.drawn_points = []
        self.labels = []
        plt.connect('button_press_event', self.onclick)

    def onclick(self, event):
        self.coordinates.append((event.xdata, event.ydata))
        p = plt.plot(event.xdata, event.ydata, 'o')
        self.drawn_points.append(p)
        plt.show()

    def clear(self):
        self.coordinates.clear()
        self.skin_regions.clear()
        self.labels.clear()
        for p in self.drawn_points:
            p.remove()
        plt.gcf().canvas.draw()

    def to_dict(self):
        seg = {}
        seg['coordinates'] = self.coordinates
        return seg

import os.path
def get_show_image(filename, band=None):
    extension = os.path.splitext(filename)[1]

    if extension == '.png':
        rgb = imageio.imread(filename)
    else:
        container = hyperread(filename)

        if band is None:
            rgb = container.rgb_image(normalize_to_minmax=True)
        else:
            rgb = container[:,:,args.band]
    return rgb

if __name__ == '__main__':
    #parse arguments
    parser = argparse.ArgumentParser(description='Select regions or point coordinates in images.')
    parser.add_argument('files', metavar='file', type=str, nargs='+',
            help='Path to hyperspectral image or PNG image')
    parser.add_argument('--output-path', metavar='output_path', type=str, default='segmentation_list.p',
            help='Output path.')
    parser.add_argument('--append', action='store_true', default=False,
            help='Whether to append data to existing file.')
    parser.add_argument('--band', type=int, default=None, help='Band to display instead of RGB.')
    parser.add_argument('--points', default=False, action='store_true', help='Select points instead of regions')

    args = parser.parse_args()

    if args.append:
        seg_list = pickle.load(open(args.output_path, 'rb'))
    else:
        seg_list = []

    for i, filename in enumerate(args.files):
        rgb = get_show_image(filename, args.band)

        #show image
        plt.imshow(rgb)
        ax = plt.gca()

        #select regions
        if args.points:
            handler = point_selection_handler()
        else:
            handler = region_selection_handler()
        plt.ioff()
        plt.show()

        seg = handler.to_dict()
        seg['filename'] = filename
        seg_list.append(seg)

    pickle.dump(seg_list, open(args.output_path, "wb"))
