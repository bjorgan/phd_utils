"""
Utilities for reading and writing hyperspectral images in ENVI format.

Asgeir Bjorgan, NTNU
"""

import re
import os
import numpy

class hyperwrite:
    """
    Class for writing hyperspectral files to disk.

    Attributes:
        base_header: Base header containing default bands and wavelengths, if specified.
        base_filename: Base filename for output files.
        file_id: File identificator for hyperspectral file.
        num_lines: Number of lines currently written to the hyperspectral file.
    """

    def __init__(self, base_filename, base_header=None):
        """
        Open hyperspectral file for writing. See write() for actual image data writing
        and close() for closing the file and writing a header.

        Parameters
        ----------
        base_filename: str
            Base filename for output file. Files will be saved as base_filename.img and base_filename.hdr.
        base_header: dict
            Has to contain 'wlens' and 'default_bands' for wavelengths and default RGB bands. Other fields will be ignored
            and inferred from the actual written data. If not specified, wavelengths will consist of numbers from 0 to num_bands,
            and default_bands will be the first three bands.
        """

        self.base_header = base_header
        self.base_filename = base_filename
        self.file_id = open(base_filename + '.img', 'wb')
        self.num_lines = 0

    def write(self, image_data):
        """
        Write image data to file and update internal book-keeping variables over the number of lines.

        Parameters
        ----------
        image_data: ndarray
            Image data in dimensions NUM_LINES x NUM_SAMPLES x NUM_BANDS or NUM_SAMPLES x NUM_BANDS.
        """

        #keep last written data for header inference
        self.last_image_data = image_data

        #update number of lines
        if len(image_data.shape) == 3:
            self.num_lines += image_data.shape[0]
        else:
            self.num_lines += 1

        #write to file in bil format (swap second to last and last axes)
        image_data.swapaxes(-2, -1).tofile(self.file_id)

    def close(self):
        """
        Close hyperspectral file and write header information to file.
        """
        if self.base_header is None:
            self.base_header = {}
            self.base_header['wlens'] = np.arange(0, self.last_image_data.shape[-1])
            self.base_header['default_bands'] = np.arange(0, 3).astype(int)

        #infer header properties from last line of data and accumulated lines
        header = {}
        header["samples"] = self.last_image_data.shape[-2]
        header["lines"] = self.num_lines
        header["bands"] = self.last_image_data.shape[-1]
        header["offset"] = 0
        header["interleave"] = 'bil'
        header["wlens"] = self.base_header['wlens'].copy()
        header["default_bands"] = self.base_header['default_bands'].copy()
        header["datatype"] = self.last_image_data.dtype

        #write header to file
        header_string = _hyperspectral_header_to_string(header)
        header_file = open(self.base_filename + '.hdr', 'w')
        header_file.write(header_string)
        header_file.close()

        #close hyperspectral file
        self.file_id.close()

from enum import Enum
class RgbNormalization(Enum):
    """
    Used for specifying how the RGB image constructed from the hyperspectral
    image bands should be normalized.
    """

    NO_NORMALIZATION = 0, #No normalization (assume already normalized to 0 and 1, e.g. reflectance calibrated image)
    MINMAX_NORMALIZATION = 1, #Normalize so that minimum value in image => 0, maximum value => 1
    STATS_NORMALIZATION = 2 #Normalize so that mean(img) - 2*stdev(img) => 0, mean(img) + 2*stdev(img) => 1

class hyperspectral:
    """
    Hyperspectral datacontainer, for combining image and metadata. Used as
    base class for hyperread, which maps a hyperspectral file using mmap.
    """

    def __init__(self, header, image=None):
        self.header = header
        self._image = image

    def image(self):
        return self._image

    def __getitem__(self, index):
        """
        Enable direct access through [] on the class instance.
        """
        return self.image()[index]

    def rgb_image(self, rgb_normalization = RgbNormalization.NO_NORMALIZATION):
        """
        Return RGB image constructed using the default RGB bands.

        Parameters
        ----------
        rgb_normalization: RgbNormalization instance, optional
            Whether the RGB image should be normalized, and how.
        Returns
        -------
        rgb_image: ndarray
            RGB image.
        """
        full_image = self.image()
        rgb_bands = self.header["default_bands"]
        return rgb_image(full_image, rgb_bands, rgb_normalization)

class hyperread(hyperspectral):
    """
    Read hyperspectral images as mmapped files.

    Attributes:
    filename    Image filename
    header      Image header (dict), containing image properties:
                   - "samples", the number of pixels across-track
                   - "lines", the number of pixels along-track
                   - "bands", the number of wavelengths
                   - "offset", the offset of the image data from the start of
                     the image file (hyperread takes care of this
                     automatically)
                   - "interleave", how the data is interleaved within the file
                     (the image() routine takes care of this automatically)
                   - "wlens", the wavelengths corresponding to each band index
                   - "default_bands", the default, defined RGB bands

    The most useful functions for you is probably going to be image() and
    rgb_image() (inherited from parent class).
    """

    def _read_header(self, filename):
        """
        Read hyperspectral header file, save to header dict in class.

        Parameters
        ----------
        filename: str
            Header filename
        """
        self.filename = filename
        filename = os.path.splitext(filename)[0] + '.hdr'
        fp = open(filename, "r")
        header = fp.read()

        self.header = {}
        self.header["samples"] = int(_extract_property(header, "samples"))
        self.header["lines"] = int(_extract_property(header, "lines"))
        self.header["bands"] = int(_extract_property(header, "bands"))
        self.header["offset"] = int(_extract_property(header, "header offset"))
        self.header["interleave"] = _extract_property(header, "interleave")
        self.header["wlens"] = _extract_array_property(header, "wavelength", float)
        self.header["default_bands"] = _extract_array_property(header, "default bands", int)

        datatype_str = _extract_property(header, "data type")

        self.header["datatype"] = _python_datatype_from_envi_datatype(datatype_str)
        fp.close()

    def image_raw(self):
        """
        Return memory map of hyperspectral image file. Will behave like a
        normal numpy matrix.

        Image is returned with dimensions according to the underlying
        interleave structure:

            * BIL: [LINES x BANDS x SAMPLES]
            * BIP: [LINES x SAMPLES x BANDS]
            * BSQ: [BANDS x LINES x SAMPLES]

        Use image()-function to return as [LINES x SAMPLES x BANDS] regardless
        of interleave type.

        Returns
        -------
        img: np.memmap/ndarray
            (Read-only) hyperspectral image matrix
        """
        if self.header["interleave"] == 'bil':
            shape = (self.header["lines"], self.header["bands"], self.header["samples"])
        elif self.header["interleave"] == 'bip':
            shape = (self.header["lines"], self.header["samples"], self.header["bands"])
        elif self.header["interleave"] == 'bsq':
            shape = (self.header["bands"], self.header["lines"], self.header["samples"])
        else:
            raise ValueError("Unrecognized interleave in header")

        fp = numpy.memmap(self.filename, offset=self.header["offset"], dtype=self.header["datatype"], mode='r', shape=shape)
        return fp

    def image(self):
        """
        Return hyperspectral image as LINES x SAMPLES x BANDS-matrix.

        Returns
        -------
        img: np.memmap/ndarray
            Hyperspectral image matrix
        """
        if self.header["interleave"] == 'bil':
            return numpy.swapaxes(self.image_raw(), 1, 2)
        elif self.header["interleave"] == 'bip':
            return self.image_raw()
        elif self.header["interleave"] == 'bsq':
            return numpy.swapaxes(self.image_raw(), 0, 2)
        else:
            raise ValueError("Unrecognized interleave in header")

    def __init__(self, filename):
        """
        Prepare image header and filename of hyperspectral image. Use
        .image*()-routines to extract the actual hyperspectral image as a numpy
        matrix, or access the data using the []-operator directly on the class
        instance.

        Parameters
        ----------
        filename: str
            Input hyperspectral image filename
        """
        self._read_header(filename)

def rgb_image(image, rgb_bands, rgb_normalization=RgbNormalization.NO_NORMALIZATION):
    """
    Construct RGB image from hyperspectral image data.

    Parameters
    ----------
    image: ndarray (LINES x SAMPLES x BANDS)
        Image data
    rgb_bands: list
        Band indices corresponding to assumed red, green and blue channels.
    rgb_normalization:
        Normalization to use, if any.
    """

    rgb = image[:,:,rgb_bands].copy()

    if rgb_normalization == RgbNormalization.MINMAX_NORMALIZATION:
        max_values = numpy.max(numpy.max(rgb, axis=0), axis=0)
        min_values = numpy.min(numpy.min(rgb, axis=0), axis=0)

    if rgb_normalization == RgbNormalization.STATS_NORMALIZATION:
        meanval = numpy.mean(numpy.mean(rgb, axis=0), axis=0)
        stdval = numpy.std(numpy.std(rgb, axis=0), axis=0)
        max_values = meanval + 2*stdval;
        min_values = meanval - 2*stdval;

    if rgb_normalization is not RgbNormalization.NO_NORMALIZATION:
        max_values = max_values[None,None,:]
        min_values = min_values[None,None,:]
        rgb = (rgb - min_values)/(max_values - min_values)

    rgb[rgb < 0] = 0
    rgb[rgb > 1] = 1
    return rgb

def _hyperspectral_header_to_string(header):
    """
    Convert hyperspectral header information to valid ENVI string, for writing to header file.

    Parameters
    ----------
    header: dict
        ENVI header information
    Returns
    -------
    header_text: str
        ENVI header string
    """

    header_text = 'ENVI\nsamples = ' + str(header['samples'])
    header_text += '\nlines = ' + str(header['lines'])
    header_text += '\nbands = ' + str(header['bands'])
    header_text += '\nheader offset = ' + str(header['offset'])
    header_text += '\nfile type = ENVI Standard'
    header_text += '\ndata type = ' + _envi_datatype_from_python_datatype(header['datatype'])
    header_text += '\ninterleave = bil'
    header_text += '\ndefault bands = {'
    for wlen in header['default_bands']:
        header_text += str(wlen) + ' '
    header_text += '}'
    header_text += '\nbyte order = 0'
    header_text += '\nwavelength = {'
    for wlen in header['wlens']:
        header_text += str(wlen) + ' '
    header_text += '}'
    return header_text

def _python_datatype_from_envi_datatype(envi_datatype):
    """
    Get appropriate Python datatype from ENVI datatype.
    """

    if envi_datatype == '12':
        datatype = numpy.dtype(numpy.uint16)
    elif envi_datatype == '4':
        datatype = numpy.dtype(numpy.float32)
    else:
        raise(ValueError("Unrecognized datatype in header file"))
    return datatype

def _envi_datatype_from_python_datatype(python_datatype):
    """
    Get appropriate ENVI datatype ID from Python datatype.
    """

    if python_datatype == numpy.dtype(numpy.uint16):
        return '12'
    elif python_datatype == numpy.dtype(numpy.float32):
        return '4'
    else:
        raise(ValueError("Unrecognized datatype in header file"))

def _extract_property(header_string, field):
    """
    Helper function for extracting property from PROPERTY_NAME = PROPERTY_VALUE style information in hyspex header files.

    Parameters
    ----------
    header_string: str
        Full header text
    field: str
        Header property
    Returns
    -------
    property_value: str
        Property value
    """
    pattern = re.compile(field + "\\s*=\\s*(.*)", re.IGNORECASE)
    match = pattern.search(header_string)
    return match.group(1)

def _extract_array_property(header_string, field, dtype=float):
    """
    Helper function for extracting array property from PROPERTY_NAME = {PROPERTY_VALUE_1, PROPERTY_VALUE_2, ...} style information in hyspex header files. Can span multiple lines.

    Parameters
    ----------
    header_string: str
        Full header text
    field: str
        Header property
    dtype: optional, datatype
        Assumed datatype
    Returns
    -------
    property_array: list
        Property array
    """

    header_string = ''.join(header_string.splitlines())

    pattern = re.compile(field + "\\s*=\\s*{\\s*([0-9|.| |,]*)}", re.IGNORECASE)
    match = pattern.search(header_string)

    if match is None:
        return None

    array_str = match.group(1)
    array = [dtype(i) for i in array_str.strip('{}').replace(',', ' ').split()]

    return array

