"""
Train random forest classifier from region lists produced using
region_selector.
"""

import argparse
import matplotlib.pyplot as plt
from hyperspectral_utils import hyperread
import matplotlib.patches as patches
import numpy as np
from matplotlib.widgets import RectangleSelector
import pickle
import random_forest

def load_regions(filename):
    region_list = pickle.load(open(filename, 'rb'))
    return region_list

def region_list_to_images_and_segmentations(filename):
    """
    Load hyperspectral data and construct segmentations
    from region lists. Chooses only the first region within the skin
    and background lists corresponding to each hsi filename.
    """

    region_list = load_regions(filename)
    filenames = np.array([file_entry['filename'] for file_entry in region_list])

    print('Warning: selecting only first region in list')
    true_coords = np.array([file_entry['skin_pos'][0] for file_entry in region_list]).astype(int)
    true_sizes = np.array([file_entry['skin_size'][0] for file_entry in region_list]).astype(int)

    false_coords = np.array([file_entry['background_pos'][0] for file_entry in region_list]).astype(int)
    false_sizes = np.array([file_entry['background_size'][0] for file_entry in region_list]).astype(int)

    image_data, segmentations = random_forest.prepare_images_and_segmentations_from_lists(filenames, true_coords, true_sizes, false_coords, false_sizes)
    return image_data, segmentations

def train_random_forest_from_region_list(filename):
    """
    Train random forest for segmentation from region list.
    """
    image_data, segmentations = region_list_to_images_and_segmentations(filename)
    rf = random_forest.hyper_rf(num_estimators=50)

    print('Training random forest')
    rf.train(image_data, segmentations)
    return rf

if __name__ == '__main__':
    #parse arguments
    parser = argparse.ArgumentParser(description='Train random forest on region list produced using region_selector.py.')
    parser.add_argument('--output-path', metavar='output_path', type=str, default='random_forest.p',
            help='Output filename for random forest file.')
    parser.add_argument('file', metavar='file', type=str, nargs=1,
            help='Path to region file.')

    args = parser.parse_args()
    output_path = args.output_path
    region_filename = args.file[0]

    rf = train_random_forest_from_region_list(region_filename)
    rf.to_file(output_path)
