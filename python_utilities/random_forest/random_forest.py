"""
Methods for training random forests to classify hyperspectral images into skin
and non-skin (or other types of segmentation tasks).

Asgeir Bjorgan, NTNU
"""

import numpy as np
import matplotlib.pyplot as plt
import sklearn.ensemble
import scipy.ndimage
from hyperspectral_utils import hyperread
import pickle

def hyper_rf_from_file(filename):
    """
    Load already trained random forest from file.

    Parameters
    ----------
    filename: str
        Filename

    Returns
    -------
    ret_rf: class hyper_rf
        Random forest
    """

    input_dict = pickle.load(open(filename, "rb"))

    ret_rf = hyper_rf()
    ret_rf.apply_preprocessing = input_dict['apply_preprocessing']
    ret_rf.apply_postprocessing = input_dict['apply_postprocessing']
    ret_rf.random_forest = input_dict['random_forest']
    return ret_rf


class hyper_rf:
    """
    Convenience wrapper around sklearn.ensemble.RandomForestClassifier for
    applying random forest to hyperspectral images for e.g. skin segmentation.

    Entry method: Use train_on_data() or train() to train the random forest on
    hyperspectral image and corresponding manually labeled segmentations, or
    load a trained random forest from file using hyper_rf_from_file() above.

    Use apply() to apply random forest to image.
    """

    def __init__(self, apply_preprocessing=True, apply_postprocessing=True, num_estimators=10):
        self.apply_preprocessing = apply_preprocessing
        self.apply_postprocessing = apply_postprocessing
        self.num_estimators=num_estimators
        self.random_forest = None
        self.region_selection = region_selection.SELECT_LARGEST

    def to_file(self, filename):
        """
        Save random forest to file.

        Parameters
        ----------
        filename: str
            Output filename
        """

        output_dict = {'random_forest': self.random_forest,
                       'apply_preprocessing': self.apply_preprocessing,
                       'apply_postprocessing': self.apply_postprocessing}
        pickle.dump(output_dict, open(filename, "wb"))

    def preprocess_data_matrix(self, image):
        """
        Sum-normalizes the spectra in the hyperspectral data matrix.

        Parameters
        ----------
        image: ndarray
            Hyperspectral image in dimensions SPATIAL_DIMENSIONS x NUM_BANDS

        Returns
        -------
        image: ndarray
            Sum-normalized image
        """
        #normalize to sum
        sum_values = np.sum(image, axis=1)
        image = image/sum_values[:, None]
        return image

    def train_on_data(self, reshaped_image, reshaped_segmentation):
        """
        Train on reshaped data arrays.

        Parameters
        ----------
        reshaped_image: ndarray, float
            Hyperspectral image array as NUM_PIXELS x NUM_BANDS matrix
        reshaped_segmentation: ndarray, boolean
            Segmentation, 1 corresponding to skin and 0 to background, length
            NUM_PIXELS
        """

        #preprocess image
        if self.apply_preprocessing:
            reshaped_image = self.preprocess_data_matrix(reshaped_image)

        #prepare and train random forest
        self.random_forest = sklearn.ensemble.RandomForestClassifier(n_estimators=self.num_estimators, n_jobs=-1)
        self.random_forest.fit(reshaped_image, reshaped_segmentation)


    def train(self, hyperspectral_images, segmentations):
        """
        Train random forest using input hyperspectral image and corresponding segmentation.

        Parameters
        ----------
        hyperspectral_images: list
            List over hyperspectral image containers (hyperread objects) or
            hyperspectral data arrays with correct spatial dimensions
        segmentations: list
            List over skin segmentations
        """
        for i, img in enumerate(hyperspectral_images):
            if isinstance(img, hyperread):
                hyperspectral_images[i] = img.image()


        k = np.shape(hyperspectral_images[0])[2]

        #reshape and combine images
        images = [ np.reshape(x, (-1, k)) for x in hyperspectral_images ]
        image = np.concatenate(images)

        segmentations = [ np.ravel(x) for x in segmentations ]
        segmentation = np.concatenate(segmentations)

        self.train_on_data(image, segmentation)


    def apply(self, image):
        """
        Apply trained random forest to hyperspectral image in order to obtain a
        skin classification.

        Parameters
        ----------
        image: {ndarray, float} or {hyperread object}
            Hyperspectral image, can either be a hyperread container or a numpy array
        Returns
        -------
        classes: ndarray, boolean
            Array over classes. 0 corresponding to background, 1 corresponding to skin
        """
        if isinstance(image, hyperread):
            image = image.image()

        reshape = False
        if len(image.shape) > 2:
            (n,m,k) = np.shape(image)
            image = np.reshape(image, (n*m, k))
            reshape = True

        if self.apply_preprocessing:
            image = self.preprocess_data_matrix(image)

        classes = self.random_forest.predict(image)

        if reshape:
            classes = np.reshape(classes, (n,m))

        if self.apply_postprocessing:
            classes = post_processing(classes, self.region_selection)

        return classes


def prepare_images_and_segmentations_from_lists(filenames, true_coords, true_sizes, false_coords, false_sizes):
    data = []
    segmentations = []
    for i, filename in enumerate(filenames):
        img = hyperread(filename).image()
        data.append(img[true_coords[i][1]:true_coords[i][1]+true_sizes[i][1], true_coords[i][0]:true_coords[i][0]+true_sizes[i][0], :])
        segmentations.append(np.ones(data[-1].shape[0:2]))

        data.append(img[false_coords[i][1]:false_coords[i][1]+false_sizes[i][1], false_coords[i][0]:false_coords[i][0]+false_sizes[i][0], :])
        segmentations.append(np.zeros(data[-1].shape[0:2]))
    return data, segmentations

from enum import Enum
class region_selection(Enum):
    """
    Types of region selection for post-processing of segmentation images.

    SELECT_ALL: Select all segmented regions.
    SELECT_CENTER: Label discontiguous regions, select the region placed in the center of the image.
    SELECT_LARGEST:  Label discontigous regions, select the largest region.
    """
    SELECT_ALL = 0
    SELECT_CENTER = 1
    SELECT_LARGEST = 2

def post_processing(segmentation, selection = region_selection.SELECT_CENTER):
    """
    Apply post-processing: Fill holes in segmentation, select specific discontiguous regions.

    Parameters
    ----------
    segmentation: ndarray
        Segmentation, values 1 or 0.
    selection: optional, region_selection
        Which region to select out of the discontiguous regions. See class
        region_selection for documentation of the alternatives.

    Returns
    -------
    filled_segmentation: ndarray
        Post-processed segmentation image.
    """

    #fill segmentation, label dis-contiguous regions
    filled_segmentation = scipy.ndimage.morphology.binary_fill_holes(segmentation)
    labeled_array, num_features = scipy.ndimage.label(filled_segmentation)

    if selection == region_selection.SELECT_CENTER:
        #find label corresponding to center of image
        c_x = int(np.floor(np.shape(filled_segmentation)[1]/2))
        c_y = int(np.floor(np.shape(filled_segmentation)[0]/2))
        label = labeled_array[c_y, c_x]
        filled_segmentation[np.invert(labeled_array == label)] = 0
    elif selection == region_selection.SELECT_LARGEST:
        #find label corresponding to the largest cluster of datapoints within the segmentation
        labels = labeled_array[filled_segmentation > 0]
        unique_labels, counts = np.unique(labels, return_counts=True)
        largest_cluster = unique_labels[np.argmax(counts)]
        filled_segmentation[np.invert(labeled_array == largest_cluster)] = 0
    return filled_segmentation

def autocrop(hyperspectral_image, segmentation):
    """
    Find bounding box of segmentation and crop hyperspectral image and mask.

    Parameters
    ----------
    hyperspectral_image: ndarray
        Hyperspectral image, as NUM_LINES x NUM_SAMPLES x NUM_BANDS
    segmentation: ndarray
        Image mask

    Returns
    -------
    hyperspectral_image: ndarray
        Cropped hyperspectral image
    segmentation: ndarray
        Cropped image mask
    """

    #find bounding box of mask image
    bbox = scipy.ndimage.measurements.find_objects(segmentation)[0]

    #crop image
    hyperspectral_image = hyperspectral_image[bbox[0].start:bbox[0].stop, bbox[1].start:bbox[1].stop, :]
    segmentation = segmentation[bbox[0].start:bbox[0].stop, bbox[1].start:bbox[1].stop]

    return hyperspectral_image, segmentation


def apply_mask_to_image(image, mask, return_as_image=True):
    """
    Apply segmentation mask to hyperspectral image.

    Parameters
    ----------
    image: ndarray, float
        Image array extracted from hyperspectral image container using image()
    mask: ndarray, boolean
        Image mask
    return_as_image: boolean
        Whether masked image should be returned as an image with masked-out
        parts set to zero, or a matrix containing only the masked spectra
    Returns
    -------
    image: ndarray, float
        Masked image array
    """
    (n,m,k) = np.shape(image)
    image = np.reshape(image, (n*m, k))
    mask = np.ravel(mask)

    if return_as_image:
        image[np.invert(mask), :] = 0
        image = np.reshape(image, (n,m,k))
    else:
        image = image[mask, :]
    return image
