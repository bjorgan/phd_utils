"""
Helper functions for generating monte carlo .mci files, for reading in .mco
files, and for running the MCML executable on a skin model as generated from
photon_transport.forward_model.
"""

import numpy as np
import os

def generate_mci_file(input_path, wlens, skin_model, num_photons, n_in=1, n_out=1):
    """
    Generate MCI file suitable for input to MCML.

    Parameters
    ----------
    input_path: str
        Basis for input and output filenames.
    wlens: str
        Wavelengths.
    skin_model: list of dicts
        Skin model, in the format as obtained from forward_model.normal_skin():
        n-length array, each element corresponding to a layer.
        Properties n, d, mua, musr, mus, g contained as dict in each element.
    num_photons:
        Number of photons to use in the simulation.
    n_in: float
        Refractive index in medium containing light source
    n_out: float
        Refractive index on other side of medium
    Returns
    -------
    full_input_path: str
        Path to .mci file
    output_filenames: list of str
        Paths to expected .mco files
    """

    input_basename = os.path.basename(input_path)

    #header
    output = []
    output.append("")
    output.append("1.0")

    #number of individual simulations
    num_wlens = len(wlens)
    output.append(str(num_wlens))
    output.append("")
    output_filenames = []

    #individual simulations
    for i in range(num_wlens):
        wlen = round(wlens[i])

        #header
        output.append("#### Wavelength " + str(wlen) + " nm")

        #output filename
        output_filename = input_basename + "_" + str(wlen) + ".mco"
        output_filenames.append(output_filename)
        output.append(output_filename + "    A")

        #number of photons
        output.append(str(num_photons))

        #discretization
        output.append("0.002 0.01")
        output.append("1 1    1")
        output.append("")

        #model
        output.append(str(len(skin_model)))
        output.append("# n    mua    mus g     d")

        #optical properties on one side of the slab
        output.append(str(round(n_in, 4)))

        #optical properties of each layer
        for layer in skin_model:
            layer_line = []
            properties = ['n', 'mua', 'mus', 'g']
            for prop in properties:
                layer_line.append(round(layer[prop][i], 4))
            layer_line.append(layer['d'])
            layer_line = [str(item) for item in layer_line]
            output.append('\t'.join(layer_line))

        #optical properties on the other side of the slab
        output.append(str(round(n_out, 4)))
        output.append("")
    output = '\n'.join(output)

    full_input_path = input_path + ".mci"
    f = open(full_input_path, "w")
    f.write(output)
    f.close()

    return full_input_path, output_filenames, output


def read_mci_file(output_filename):
    """
    Read MCML output from .mco file.
    """
    with open(output_filename, "r") as output_file:
        data = np.array(output_file.readlines())

        #diffuse reflectance
        R_d = float(data[np.char.find(data, 'Diffuse') > 0][0].split()[0])
        T_d = float(data[np.char.find(data, 'Transmittance') > 0][0].split()[0])

        #absorption as a function of depth
        A_z_start_ind = np.arange(0, len(data))[np.char.find(data, 'A_z #A') >= 0][0]+1
        A_z_end_ind = np.arange(A_z_start_ind, len(data))[np.char.equal(data[A_z_start_ind:], '\n') > 0][0]
        A_z = data[A_z_start_ind:A_z_end_ind].astype(float)

        #obtain layer properties
        num_layers = int(data[17].split()[0])
        layer_properties = []
        for i in np.arange(0, num_layers):
            layer_properties.append(np.array(data[20+i].split())[:5].astype(float))
        layer_properties = np.array(layer_properties)

        #layer depths for the cells
        dz = float(data[14].split()[0])
        z = dz*np.arange(0, len(A_z))

        #get mua in each cell
        mua = np.zeros(len(z))
        for i in range(num_layers):
            layer = num_layers - 1 - i
            mua[z < layer_properties[layer, 4]] = layer_properties[layer, 1]

        #calculate fluence rate
        fluence_rate = A_z/mua

        return {'T': T_d, 'R_d': R_d, 'A_z': A_z, 'mua': mua, 'z': z, 'fluence': fluence_rate, 'layer_props': layer_properties}

def split_skin_model(wlens, skin_model, parts):
    """
    Split a wavelength array into a given number of
    parts, and split the input skin_model into 
    several skin_models correspondingly.
    """

    #split wavelength list in a given number of parts
    wlens_collection = []
    skin_model_collection = []
    part_length = int(len(wlens)/parts)
    for i in range(parts):
        start_ind = i*part_length
        end_ind = (i+1)*part_length

        #construct new skin model for the wavelength subset
        skin_model_part = []
        for layer in skin_model:
            new_layer = {}
            for item in list(layer):
                if hasattr(layer[item], '__len__'):
                    new_layer[item] = layer[item][start_ind:end_ind]
                else:
                    new_layer[item] = layer[item]
            skin_model_part.append(new_layer)

        #add to larger arrays
        skin_model_collection.append(skin_model_part)
        wlens_collection.append(wlens[start_ind:end_ind])
    return wlens_collection, skin_model_collection

#mcml executable
MCML_EXECUTABLE_NAME = 'pymcml_mcml_executable'

#gpumcml executable
GPUMCML_EXECUTABLE_NAME = 'pymcml_gpumcml_executable'

from tqdm import tqdm
import time
import tempfile
import subprocess
from glob import glob
def _run_monte_carlo_routine(cmd_base, wlens, skin_model, n_jobs=1, num_photons=1000, n_in=1, n_out=1):
    """
    Generate .mci-files and run Monte Carlo routine. Do not use directly, use
    rather run_mcml() or run_gpumcml(), which calls this function. See their
    respective command line arguments for explanation.
    """

    #put input files in tempdir, run mcml from tempdir so that output files are in tempdir
    with tempfile.TemporaryDirectory() as tempdir:
        #split into several .mci files
        input_files = []
        output_files = []
        wlens_collection, skin_models = split_skin_model(wlens, skin_model, n_jobs)
        n_jobs = len(wlens_collection)

        for i in range(n_jobs):
            input_file_job, output_files_job, _ = generate_mci_file(tempdir + '/input' + str(i), wlens_collection[i], skin_models[i], num_photons, n_in=n_in, n_out=n_out)
            output_files.append(output_files_job)
            input_files.append(input_file_job)
        output_files = np.concatenate(output_files)

        #run MCML
        procs = []
        for i in range(n_jobs):
            #batch simulation launched in forked process, run from tempdir so that output files are put there
            cmd = cmd_base.copy()
            cmd.append(input_files[i])
            procs.append(subprocess.Popen(cmd, cwd=tempdir, stdout=open('/dev/null', 'w')))

        #check which processes are finished, update a progress bar
        if n_jobs > 1:
            pbar = tqdm(total=n_jobs)

        not_finished = procs
        finished = []
        while True:
            for i, proc in enumerate(not_finished):
                if proc.poll() is not None:
                    finished.append(not_finished.pop(i))
                    if n_jobs > 1:
                        pbar.update(1)
                    break
            time.sleep(0.1)

            if len(not_finished) == 0:
                break

        #get output
        output = []
        for output_file in np.sort(output_files):
            output.append(read_mci_file(tempdir + "/" + output_file))

        R_d = [o['R_d'] for o in output]
        T = [o['T'] for o in output]

        if n_jobs > 1:
            pbar.close()

        return {'R_d': R_d, 'T': T}

def run_mcml(wlens, skin_model, n_jobs=1, num_photons=1000, n_in=1, n_out=1):
    """
    Run MCML on supplied skin model.

    Parameters
    ----------
    wlens: array
        Wavelengths
    skin_model: array of dicts
        Skin model, in the format as obtained from forward_model.normal_skin():
        n-length array, each element corresponding to a layer.
        Properties n, d, mua, musr, mus, g contained as dict in each element.
    n_jobs: optional, int
        Number of parallel CPU tasks.
    num_photons: optional, int.
        Number of photons to use in the simulation.
    n_in: float
        Refractive index in medium containing light source
    n_out: float
        Refractive index on other side of medium
    use_gpumcml: optional, boolean
        Whether to use GPU-MCML instead of MCML. Forces the number of jobs to 1.
    saferprimes_path: str, optional
        Path to safeprimes data file. Needs to be supplied.  Will be installed
        at INSTALL_PREFIX/local/gpumcml/safeprimes_base32.txt.  Might rather
        have been found automatically, but a bit of an hassle.

    Returns
    -------
    output_dict: dict
        Contains R_d as the diffuse reflectance, T as the transmittance.
    """

    cmd_base = [MCML_EXECUTABLE_NAME]
    return _run_monte_carlo_routine(cmd_base, wlens, skin_model, n_jobs, num_photons, n_in, n_out)

def run_gpumcml(wlens, skin_model, safeprimes_path, num_photons=1000, n_in=1, n_out=1, detect_absorption=False):
    """
    Run GPU-MCML on supplied skin model.

    Parameters
    ----------
    wlens: array
        Wavelengths
    skin_model: array of dicts
        Skin model, in the format as obtained from forward_model.normal_skin():
        n-length array, each element corresponding to a layer.
        Properties n, d, mua, musr, mus, g contained as dict in each element.
    saferprimes_path: str
        Path to safeprimes data file. Needs to be supplied.  Will be installed
        at INSTALL_PREFIX/local/gpumcml/safeprimes_base32.txt.  Might rather
        have been found automatically, but a bit of an hassle.
    num_photons: optional, int.
        Number of photons to use in the simulation.
    n_in: float
        Refractive index in medium containing light source
    n_out: float
        Refractive index on other side of medium
    detect_absorption: boolean, optional
        Whether to record absorption. Speeds up simulation a bit if set to False.

    Returns
    -------
    output_dict: dict
        Contains R_d as the diffuse reflectance, T as the transmittance.
    """

    cmd_base = [GPUMCML_EXECUTABLE_NAME, '-p' + safeprimes_path]
    if not detect_absorption:
        cmd_base.append('-A')
    n_jobs = 1

    return _run_monte_carlo_routine(cmd_base, wlens, skin_model, n_jobs, num_photons, n_in, n_out)
