"""
Routines for applying the MNF transform to an image or multiple images.
"""

import numpy as np
import scipy.linalg
from tqdm import tqdm

class mnf:
    """
    Convenience class for obtaining the MNF transform from an input image.

    The same functionality can also be obtained by running the individual
    functions defined in this module below.

    This class has functionality for orthonormalizing the transform. Some notes on this:

    If the first n bands in the inverse transform are
    orthonormalized (and forward transform
    modified accordingly), we have some nice properties.

    Preservation of euclidean distance
    ----------------------------------

    Assume that denoised spectrum y (B x 1) would be

    y = A z + mu,

    where A is the first n bands of the inverse MNF transform (B x n), z is the n first bands of the transformed image pixel (n x 1), mu is the band mean (B x 1). First n bands of inverse are orthogonalized, so that A^T A = I.

    Then the Euclidean distance squared between y_1 and y_2 can be expressed as

    ||(y_1 - y_2)||^2
    = (y_1 - y_2)^T (y_1 - y_2)
    = (A z_1 + mu - A z_2 - mu)^T (A z_1 + mu - A z_2 - mu)
    = (z_1 - z_2)^T A^T A (z_1 - z_2)
    = (z_1 - z_2)^T (z_1 - z_2)
    = ||(z_1 - z_2) ||^t

    This means that distances between denoised spectra
    are the same as the distances between the spectra
    expressed in the MNF transformed space. Euclidean
    distance can therefore be calculated directly.

    Dot product
    -----------
    Dot product (and length of a vector) are a bit more
    of an hassle.

    y_1 \\dot y_2
    = y_1^T y_2
    = (A z_1 + mu)^T (A z_2 + mu)
    = z_1^T A^T A z_2 + z_1^T A^T mu + mu^T A z_2 + mu^T mu
    = z_1^T z_2 + mu^T A (z_1 + z_2) + mu^T mu

    Dot product between y_1 and y_2 is expressed as dot product
    between z_1 and z_2, but with some extra stuff. Here, mu^T A and ||mu||^2
    can be precomputed. See dot_product() for an
    implementation of this.
    """

    def __init__(self, image, num_bands_in_inverse, orthonormalize=True, verbose=True):
        """
        Create MNF object.

        Parameters
        ----------
        image: ndarray, float
            Input image of dimensions lines x samples x bands.
        num_bands_in_inverse: int
            Number of bands which would be used in the inverse
        orthonormalize: boolean, optional
            Whether to orthonormalize the assumed noise-free components
        verbose: boolean, optional
            Enable verbose output (progress bars, ...)
        """

        self.verbose = verbose
        self.num_bands_in_inverse = num_bands_in_inverse

        #statistics
        self.cov = estimate_image_covariance(image, verbose)
        self.noise_cov = estimate_noise_covariance(image, verbose)
        self.orthonormalized = orthonormalize

        #calculate transforms
        if orthonormalize:
            self.orthonormalize_bands = num_bands_in_inverse
        else:
            self.orthonormalize_bands = None

        self.recalculate_transforms()

    def recalculate_transforms(self):
        """
        Update forward and inverse transform with current covariances.
        """
        self.forward, self.inverse, self.eigenvalues = mnf_transformation_matrix(self.cov, self.noise_cov, orthonormalize_bands=self.orthonormalize_bands)

    def transform(self, image):
        """
        Calculate forward MNF transform of input image.

        Parameters
        ----------
        image: ndarray, float
            Image of dimensions lines x samples x bands.
        Returns
        -------
        transformed:
            Transformed image
        """
        transformed, _ = transform_forward(self.forward, image, keep_num_bands=self.num_bands_in_inverse, verbose=self.verbose, mean_spectrum=self.cov.get_mean())
        return transformed

    def add_covariances(self, image_cov, noise_cov):
        """
        Update covariances with new covariance information and recalculate
        forward and inverse transforms.

        Parameters
        ----------
        image_covariance: class covariance instance
            Image covariance, obtained from estimate_image_covariance().
        noise_covariance: class covariance instance
            Noise covariance, obtained from estimate_noise_covariance().
        """
        self.cov.update_with_covariance(image_cov)
        self.noise_cov.update_with_covariance(noise_cov)
        self.recalculate_transforms()

    def dot_product(self, transformed_vec_1, transformed_vecs_2):
        """
        Calculate dot product between denoised spectra
        from spectra in MNF transformed coordinates, subset
        to the number of bands that would have been
        used in inverse.

        See docstring of this class for explanation
        of how the dot product can be calculated in this way.

        Parameters
        ----------
        transformed_vec_1:
            MNF transformed spectrum, n x 1 (ordinary 1D array)
        transformed_vecs_2:
            MNF transformed spectra, n x (arbitrary number).

        Returns
        -------
        dot_product:
            Dot product as if dot product was calculated between
            the denoised spectra in the original wavelength space.
        """
        if not self.orthonormalized:
            raise("Can't calculate dot products without orthornormalized inverse transform")

        row_vec = transformed_vec_1[:, np.newaxis]
        mean_vec = self.cov.get_mean()[:, np.newaxis]

        if len(transformed_vecs_2.shape) == 1:
            transformed_vecs_2 = transformed_vecs_2[:, np.newaxis]

        dot_product = row_vec.T @ transformed_vecs_2 + mean_vec.T @ self.inverse[:, 0:self.num_bands_in_inverse] @ (row_vec + transformed_vecs_2) + mean_vec.T @ mean_vec
        return np.squeeze(dot_product)

    def vector_length(self, transformed):
        """
        Calculate euclidean length of input transformed image
        as if they were in denoised spectra land. See dot_product() docstring
        for more details.

        Parameters
        ----------
        transformed:
            Transformed image, in dimensions lines x samples x bands.
        """
        reshaped = transformed.reshape((-1, transformed.shape[2])).T
        if len(reshaped.shape) == 1:
            reshaped = reshaped[:, np.newaxis]
        mean_vec = self.cov.get_mean()[:, np.newaxis]

        #should ideally have reused dot_product function, but hassle
        #length_squared = reshaped.T @ reshaped + mean_vec.T @ self.inverse[:, 0:self.num_bands_in_inverse] @ (2*reshaped) + mean_vec.T @ mean_vec
        length_squared = np.sum(reshaped**2, axis=0) + (mean_vec.T @ self.inverse[:, 0:self.num_bands_in_inverse]) @ (2*reshaped) + mean_vec.T @ mean_vec

        length_squared = length_squared.reshape(transformed.shape[0:2])

        return np.squeeze(np.sqrt(length_squared))

    def denoise(self, image):
        """
        Calculate denoised spectra using forward and inverse MNF transforms.

        Parameters
        ----------
        image: ndarray, float
            Image of dimensions lines x samples x bands.
        Returns
        -------
        denoised:
            Denoised image
        """
        return denoise(self.forward, self.inverse, self.num_bands_in_inverse, image)

def orthonormalize(forward, inverse, num_bands_in_inverse):
    """
    Orthornormalize the first n components of the inverse MNF transform,
    and modify the forward transform accordingly, so that
    the first n bands of the MNF transform represent coordinates
    in an orthonormal basis.

    See docstring of orthonormalize_bands argument in mnf_transformation_matrix()
    for motivation and more information.

    Usually no need to call this function explicitly since
    mnf_transformation_matrix() has the orthonormalize_bands argument to do the
    orthonormalization implicitly.

    Parameters
    ----------
    forward:
        MNF forward transformation matrix as obtained from mnf_transformation_matrix().
    inverse:
        MNF inverse transformation matrix as obtained from mnf_transformation_matrix().

    num_bands_in_inverse:
        Number of bands to orthonormalize from start of the transform.

    Returns
    -------
    new_forward:
        New forward transformation matrix, modified to let first n components
        be orthonormal with respect to each other
    new_inverse:
        Corresponding inverse transformation matrix
    """

    #do QR decomposition of the first n bands we would use in the inverse transform
    #(Q has orthonormal column vectors)
    Q, R = np.linalg.qr(inverse[:, :num_bands_in_inverse])

    #extend Q and R so that Q * R = inverse, and not only the first part that was decomposed
    Q_extended = inverse.copy()
    Q_extended[:, :num_bands_in_inverse] = Q

    R_extended = np.eye(len(Q_extended))
    R_extended[:num_bands_in_inverse, :num_bands_in_inverse] = R

    #create new forward and inverse matrices
    new_forward = R_extended @ forward

    #(new_forward)^-1 is now forward^-1 R^-1 = inverse R^-1 = Q R R^-1 = Q,
    #so that Q is the inverse of new_forward
    new_inverse = Q_extended

    return new_forward, new_inverse


def mnf_transformation_matrix(image_covariance, noise_covariance, orthonormalize_bands=None):
    """
    Calculate MNF transformation matrices.

    Parameters
    ----------
    image_covariance: class covariance instance
        Image covariance, obtained from estimate_image_covariance().
    noise_covariance: class covariance instance
        Noise covariance, obtained from estimate_noise_covariance().
    orthonormalize_bands: int, optional
        Whether to orthononalize the forward transform. If set to None or 0:
        Don't orthonormalize. If set to an integer n > 0: Orthonormalize the
        first n bands in the transform. n would be equivalent with the number
        of bands to use in the inverse transform for denoising.

        This means that the first (noise-free) components are orthonormalized,
        while the rest of the (noisy) components are left as they are. This
        would be desired if mathematical operations on MNF-transformed
        coordinates are to yield equivalent results with mathematical
        operations on original coordinates, since this is an orthogonal basis.
        Noise bands are not orthornormalized as they would otherwise just be
        skipped.

        The first n bands of the the inverse transform is orthonormalized, and
        the forward transform is modified accordingly so that the resulting
        inverse transform is its inverse.

        After orthonormalization, SNR ordering within the orthonormalized bands
        is no longer correct since the bands will be mixed, but the assumption
        is that these bands would be treated equally in any case.

    Returns
    -------
    forward_transform: ndarray, float
        MNF forward transform
    inverse_transform: ndarray, float
        MNF inverse transform.
    eigenvalues: ndarray, float
        Eigenvalues.
    """

    #solve generalized symmetric eigenvalue problem: noise_covariance * eigvec = lambda * image_covariance * eigvec
    eigenvalues, eigenvectors = scipy.linalg.eigh(noise_covariance.cov(), image_covariance.cov(), type=1)

    forward = eigenvectors.T
    inverse = np.linalg.inv(forward)

    if orthonormalize_bands is not None:
        forward, inverse = orthonormalize(forward, inverse, orthonormalize_bands)

    #inverse transform is the inverse of the forward transform
    return forward, inverse, eigenvalues

def estimate_noise_covariance(image, verbose=False):
    """
    Estimate noise covariance of input image. Assumes spatial correlation,
    so that noise is estimate as the pixel difference.

    Parameters
    ----------
    image: ndarray
        Hyperspectral image

    Returns
    -------
    noise_covariance: class covariance instance
        Noise covariance matrix
    """

    if verbose:
        pbar = tqdm(total=image.shape[0], desc='Estimate noise covariance')

    cov = covariance(image.shape[2])
    for i in np.arange(0, image.shape[0]):
        line = image[i, :, :]
        noise_estimate = line[0:-1, :].astype(float) - line[1:, :].astype(float)
        cov.update_with_image(noise_estimate)

        if verbose:
            pbar.update(1)

    if verbose:
        pbar.close()
    return cov

class covariance:
    """
    Covariance structure for adding covariances in a numerically stable way.

    Add data to covariance using update_with_image(), or combine the current
    covariance with another covariance structure using
    update_with_covariance().

    We could just have obtained the covariances directly with np.cov(), but
    seems to be a bit memory intensitive with large arrays, and the ability to
    combine covariances from multiple images is also desired for e.g. k means
    clustering (?).
    """

    def __init__(self, shape):
        self.curr_cov = np.zeros((shape, shape))
        self.curr_mean = np.zeros(shape)
        self.N_curr = 0

    def covariance_from_data(image_matrix):
        """
        Create covariance object from data.

        Parameters
        ----------
        image_matrix: ndarray, float
            Hyperspectral image matrix, of dimensions SAMPLES x BANDS.

        Returns
        -------
        cov: covariance object
            Covariance object with statistics calculated from the input image matrix
        """
        cov = covariance(image_matrix.shape[1])
        cov.N_curr = image_matrix.shape[0]
        cov.curr_cov = np.cov(image_matrix, rowvar=False, bias=True)*cov.N_curr
        cov.curr_mean = np.mean(image_matrix, axis=0)
        return cov

    def update_with_image(self, image_matrix):
        """
        Update covariance with input data.

        Parameters
        ----------
        image_matrix: ndarray, float
            Hyperspectral image matrix, of dimensions SAMPLES x BANDS.
        """
        image_matrix_stats = covariance.covariance_from_data(image_matrix)
        self.update_with_covariance(image_matrix_stats)

    def update_with_covariance(self, cov):
        """
        Update covariance with covariance information from another structure,
        so that the total result is as if the covariance was calculated from
        all of these data.

        Parameters
        ----------
        cov: covaraince object
            Covariance object calculated from some data.
        """
        N_prev = self.N_curr
        prev_mean = self.curr_mean.copy()

        added_mean = cov.curr_mean
        added_cov = cov.curr_cov
        N_added = cov.N_curr

        self.N_curr = N_added + N_prev

        mean_diff = (added_mean - prev_mean)[:, np.newaxis]
        correction = mean_diff @ mean_diff.T * N_added*N_prev/self.N_curr

        self.curr_cov = self.curr_cov + added_cov + correction
        self.curr_mean = prev_mean + N_added/self.N_curr*(added_mean - prev_mean)

    def cov(self):
        """
        Return numerical covariance matrix.

        Returns
        -------
        cov: ndarray, float
            Symmetric matrix containing the covariance values.
        """
        return self.curr_cov/self.N_curr

    def get_mean(self):
        """
        Return mean.

        Returns
        -------
        mean: ndarray float
            Array containing current mean values.
        """
        return self.curr_mean.copy()

def estimate_image_covariance(image, verbose=False):
    """
    Estimate covariance of input image.

    Parameters
    ----------
    image: ndarray
        Hyperspectral image

    Returns
    -------
    covariance: class covariance instance
        Covariance matrix
    """

    cov = covariance(image.shape[2])

    if verbose:
        pbar = tqdm(total=image.shape[0], desc='Estimate image covariance')

    for i in np.arange(0, image.shape[0]):
        cov.update_with_image(image[i, :, :])

        if verbose:
            pbar.update(1)
    if verbose:
        pbar.close()

    return cov


def denoise(forward, inverse, num_bands_in_inverse, image, mean_spectrum=None):
    """
    Denoise an image using the MNF transform.

    Parameters
    ----------
    forward: ndarray, float
        MNF forward transform, as obtained from mnf_transformation_matrix().
    inverse: ndarray, float
        MNF inverse transform, as obtained from mnf_transformation_matrix().
    num_bands_in_inverse: int
        Number of bands in inverse.
    image: ndarray, float
        Hyperspectral image

    Returns
    -------
    denoised: ndarray, float
        Denoised hyperspectral image
    """
    forward_transform, mean_spectrum = transform_forward(forward, image, mean_spectrum=mean_spectrum)
    inverse_transform = transform_inverse(inverse, forward_transform, mean_spectrum, num_bands_in_inverse)
    return inverse_transform

def transform_forward(forward, image, keep_num_bands=None, verbose=False, mean_spectrum=None):
    """
    Calculate MNF forward transform of input image.

    Parameters
    ----------
    forward: ndarray, float
        MNF forward transformation matrix, as obtained from mnf_transformation_matrix().
    image: ndarray, float
        Hyperspectral image
    keep_num_bands: int, optional
        Number of bands to keep in the transformed image. None means all
        bands are returned.
    verbose: boolean, optional
        Print progress bar
    mean_spectrum: ndarray, float, optional
        Mean to subtract from image before applying transform. Calculated from
        the image if set to None.

    Returns
    -------
    forward_transform: ndarray, float
        MNF forward transform of input image.
    mean_spectrum: ndarray, float
        Band means of input image.
    """
    if mean_spectrum is None:
        mean_spectrum = np.mean(image, axis=(0,1))

    if verbose:
        pbar = tqdm(total=image.shape[0])

    if keep_num_bands is None:
        keep_num_bands = image.shape[2]

    ret_transform = np.zeros((image.shape[0], image.shape[1], keep_num_bands))

    for line in np.arange(image.shape[0]):
        line_image = image[line, :, :]
        line_image = line_image - mean_spectrum
        ret_transform[line, :, :keep_num_bands] = (forward @ line_image.T).T[:, :keep_num_bands]

        if verbose:
            pbar.update(1)
    if verbose:
        pbar.close()

    return ret_transform, mean_spectrum

def transform_inverse(inverse, forward_transform, mean_spectrum, num_bands_in_inverse):
    """
    Calculate MNF inverse transform of MNF forward transformed image.

    Parameters
    ----------
    inverse: ndarray, float
        MNF inverse transform, as obtained from mnf_transformation_matrix().
    forward_transform: ndarray, float
        Forward transform image, as obtained from transform_forward().
    mean_spectrum: ndarray, float
        Band means of input image, as obtained from transform_forward().
    num_bands_in_inverse: int
        Number of bands in inverse.

    Returns
    -------
    inverse_transform: ndarray, float
        Inverse MNF transform of forward transformed image, with mean added. In practice, this is a denoised image.
    """
    reshaped = np.reshape(forward_transform, (-1, forward_transform.shape[2]))
    R = np.zeros(inverse.shape)
    np.fill_diagonal(R[0:num_bands_in_inverse, 0:num_bands_in_inverse], 1)

    inverse_transform = inverse @ R @ reshaped.T
    inverse_transform = inverse_transform.T + mean_spectrum
    return np.reshape(inverse_transform, forward_transform.shape)
