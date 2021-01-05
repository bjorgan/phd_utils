#!/usr/bin/env python

"""
Utilities for parsing OceanOptics ProcSpec files.

read_procspec: Reads a single ProcSpec file.
read_procspec_in_folder: Read all files in a folder.
read_procspec_files: Reads all ProcSpec files corresponding to a file pattern.
"""

import re
import os
import numpy
from dateutil import parser as dateutil_parser

def read_procspec(filename, start_wlen=400, end_wlen=850):
    """
    Read OceanOptics .ProcSpec file (tab delimited, ASCII format).

    Example:

    `reflectance, wavelengths, data_header = read_procspec('test.ProcSpec')`

    Parameters
    ----------
    filename: str
        Filename
    start_wlen: float
        Start wavelength (ignore data values below this wavelength)
    end_wlen: float
        End wavelength (ignore data values above this wavelength)

    Returns
    -------
    spectral_data: ndarray
        Spectral data with values from 0 to 1.0
    wavelengths: ndarray
        Wavelength values corresponding to each spectral datapoint
    header: dict
        Selected header contents
    """
    fp = open(filename, "r")
    contents = fp.read()
    date_string = re.search("Date: (.*)", contents).group(1)
    header = {"date": dateutil_parser.parse(date_string)}

    #remove header and footer
    contents = re.compile("^.*>{5}Begin.*<{5}(.*)>{5}.*End", re.DOTALL).search(contents).group(1)

    #replace comma by dot
    contents = contents.replace(",", ".")

    #replace - by nothing as it fucks up things, don't matter anyway
    contents = contents.replace("-", "")

    contents = numpy.fromstring(contents, sep=" ")
    contents = numpy.reshape(contents, [int(len(contents)/2), 2])

    #limit wavelength and reflectance array to input wavelength limits
    wavelengths = contents[:,0]
    reflectance = contents[:,1]/100.0
    wavelength_limit = (wavelengths > start_wlen) & (wavelengths < end_wlen)
    wavelengths = wavelengths[wavelength_limit]
    reflectance = reflectance[wavelength_limit]

    return reflectance, wavelengths, header

import pandas as pd

def read_procspec_in_folder(folder, extension='.ProcSpec'):
    """
    Read all OceanOptics ProcSpec files in specified folder path. It is assumed that all files
    have exactly the same wavelengths.

    Example: data, wavelengths = read_procspec_in_folder('/tmp/tjafs', '.transmission')

    Parameters
    ----------
    folder: str
        Folder path
    extenstion: str
        File extension

    Returns
    -------
    ret_frame: pandas.DataFrame
        Pandas dataframe containing the spectral data, a datetime as obtained from the header and the filename.
    wavelengths: ndarray
        Wavelengths corresponding to each spectral datapoint
    """
    ret_frame = pd.DataFrame()
    for name in os.listdir(folder):
        if name.endswith(extension):
            data, wavelengths, header = read_procspec(os.path.join(folder, name))
            ret_frame = ret_frame.append(pd.Series({'spectrum': data, 'datetime': header['date'], 'filename': name}), ignore_index=True)

    return ret_frame, wavelengths

from glob import glob

def read_procspec_files(file_pattern):
    """
    Read all OceanOptics ProcSpec files using the specified file pattern. It is
    assumed that all files have exactly the same wavelengths.

    Example: data, wavelengths = read_procspec_files('/tmp/tjafs/*.transmission')

    Parameters
    ----------
    file_pattern: str
        File pattern

    Returns
    -------
    ret_frame: pandas.DataFrame
        Pandas dataframe containing the spectral data, a datetime as obtained from the header and the filename.
    wavelengths: ndarray
        Wavelengths corresponding to each spectral datapoint
    """
    ret_frame = pd.DataFrame()
    for filename in glob(file_pattern):
        data, wavelengths, header = read_procspec(filename)
        ret_frame = ret_frame.append(pd.Series({'spectrum': data, 'datetime': header['date'], 'filename': filename}), ignore_index=True)

    return ret_frame, wavelengths
