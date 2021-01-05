"""
Functions for calibration of specim hyperspectral images.  Used by executable
specim_hyper_calibrate, which writes to file. For use in other python modules,
see function calibrate_to_memory().

Asgeir Bjorgan, NTNU
"""
from hyperspectral_utils import hyperread, hyperwrite, hyperspectral
import os
import numpy as np
from tqdm import tqdm


def calibrate_image(image, white_ref, dark_ref):
    """
    Calibrate hyperspectral image according to white and dark reference standards.

    If white_ref - dark_ref is zero, the resulting calibrated value is set to 1.0.

    Parameters
    ----------
    image:
        Hyperspectral image array ([num_lines x] num_samples x num_bands)
    white_ref:
        White reference array (num_samples x num_bands)
    dark_ref:
        Dark reference array (num_samples x num_bands)
    """

    calibrated = (image - dark_ref)/(white_ref - dark_ref)
    calibrated[np.invert(np.isfinite(calibrated))] = 1.0
    return calibrated

def read_specim_files(filepath):
    """
    Read hyperspectral files produced using a SpecIm camera. Assumes the
    following files, assuming that filepath consists of dir/basename:
     - dir/basename: Path to hyperspectral file
     - dir/DARKREF_basename: Path to hyperspectral file containing dark
       reference
     - dir/WHITEREF_basename: Path to hyperspectral file containing white
       reference standard

    Parameters
    ----------
    filepath: str
        Path to file containing actual image. Dark/white reference images are
        inferred from this path according to rules above.

    Returns
    -------
    image: hyperread object
        Hyperspectral image
    white_ref: hyperread object
        White reference image
    dark_ref: hyperread object
        Dark reference image
    """

    base_filename = os.path.basename(filepath)
    dirname = os.path.dirname(filepath)

    image = hyperread(filepath)
    dark_ref = hyperread(os.path.join(dirname, 'DARKREF_' + base_filename))
    white_ref = hyperread(os.path.join(dirname, 'WHITEREF_' + base_filename))

    return image, white_ref, dark_ref

def prepare_white_and_dark_reference(white_ref_arr, dark_ref_arr):
    #take mean over y axis, but keep spatial across-track dimension
    white_ref = np.mean(white_ref_arr, axis=0).astype(np.float32)
    dark_ref = np.mean(dark_ref_arr, axis=0).astype(np.float32)

    return white_ref, dark_ref

def calibrate_and_write(filename, output_filename):
    """
    Read hyperspectral image, calibrate and write directly to file.

    Parameters
    ----------
    filename: str
        Path to hyperspectral image. Dark and white reference files are
        inferred from the filepath.
    output_filename: str
        Output filename for output filename.
    """

    image, white_ref, dark_ref = read_specim_files(filename)
    white_ref, dark_ref = prepare_white_and_dark_reference(white_ref.image(), dark_ref.image())

    #calibrate and write to file
    out = hyperwrite(output_filename, image.header)
    for line in tqdm(image.image()):
        out.write(calibrate_image(line, white_ref, dark_ref))
    out.close()

def calibrate_to_memory(filename, verbose=True, bands=None):
    """
    Read specim hyperspectral image files and return calibrated reflectance image.
    In opposition to calibrate_and_write(), which writes directly to file.

    Parameters
    ----------
    filename: str
        Path to hyperspectral image. Dark and white reference files are
        inferred from the filepath.
    verbose: boolean, optional
        Whether to write messages and progress bar to stdout
    bands: ndarray, int, optional
        If set, returns only the specific bands in calibrated image

    Returns
    -------
    ret_image: class hyperspectral object
        Reflectance image as class hypespectral object
    """

    #read and calibrate image
    image_container, white_ref, dark_ref = read_specim_files(filename)
    white_ref = white_ref.image()
    dark_ref = dark_ref.image()
    image = image_container.image()

    if bands is not None:
        image = image[:,:,bands]
        white_ref = white_ref[:,:,bands]
        dark_ref = dark_ref[:,:,bands]

    white_ref, dark_ref = prepare_white_and_dark_reference(white_ref, dark_ref)

    calibrated = np.zeros(image.shape, dtype=np.float32)

    if verbose:
        for i, line in enumerate(tqdm(image, desc='Calibrate')):
            calibrated[i,:,:] = calibrate_image(line, white_ref, dark_ref)
    else:
        calibrated = calibrate_image(image, white_ref, dark_ref)

    #modify datatype of header
    header = dict(image_container.header)
    header['datatype'] = calibrated.dtype

    #return as encapsulated image
    ret_image = hyperspectral(header, calibrated)
    return ret_image
