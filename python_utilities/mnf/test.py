"""
Tests for mnf.py. Sanity test of forward and inverse transform of a dummy
image with added noise.
"""

import unittest
import mnf
import numpy as np

class test_mnf(unittest.TestCase):
    def setUp(self):
        #create dummy image
        self.image = np.zeros((139, 145, 21))

        #spectra: linearly increasing, cosine, sine, linearly decreasing
        dummy_signal_1 = np.arange(0, self.image.shape[2])
        dummy_signal_2 = np.cos(dummy_signal_1*0.5)
        dummy_signal_3 = np.sin(dummy_signal_1*0.3)

        dummy_signal_1 = dummy_signal_1/self.image.shape[2]
        dummy_signal_4 = 1-dummy_signal_1

        #divide image in four squares, assign each of the spectra to the four boxes
        square_y = int(self.image.shape[0]/2.0)
        square_x = int(self.image.shape[1]/2.0)

        self.image[0:square_y, 0:square_x, :] = dummy_signal_1
        self.image[0:square_y, square_x:, :] = dummy_signal_2
        self.image[square_y:, 0:square_x, :] = dummy_signal_3
        self.image[square_y:, square_x:, :] = dummy_signal_4

        #add noise to image
        self.noise_stdev = 0.1
        self.noise = np.random.normal(0, self.noise_stdev, self.image.shape)
        self.noisy_image = self.image + self.noise

        #calculate image and noise covariances
        self.image_covariance = mnf.estimate_image_covariance(self.noisy_image)
        self.noise_covariance = mnf.estimate_noise_covariance(self.noisy_image)

        #calculate mnf transform
        self.forward, self.inverse, self.eigenvalues = mnf.mnf_transformation_matrix(self.image_covariance, self.noise_covariance)

    def test_forward_and_full_inverse_yields_identity(self):
        full_transform = self.inverse @ self.forward
        identity = np.eye(full_transform.shape[0])

        for i, value in enumerate(np.ravel(full_transform)):
            self.assertAlmostEqual(value, np.ravel(identity)[i])

    def test_denoised_image_is_close_to_original_noisefree_image(self):
        denoised = mnf.denoise(self.forward, self.inverse, 3, self.noisy_image)
        np.testing.assert_allclose(self.image, denoised, atol=0.5, rtol=0)

    def test_denoised_image_with_almost_full_inverse_yields_original_noisy_image(self):
        denoised = mnf.denoise(self.forward, self.inverse, self.noisy_image.shape[2]-1, self.noisy_image)
        np.testing.assert_allclose(self.noisy_image, denoised, atol=0.5, rtol=0)

    def test_denoised_image_with_full_inverse_yields_original_noisy_image(self):
        denoised = mnf.denoise(self.forward, self.inverse, self.noisy_image.shape[2], self.noisy_image)

        np.testing.assert_array_almost_equal(self.noisy_image, denoised)

    def test_image_covariance_yields_the_direct_covariance(self):
        reshaped = np.reshape(self.noisy_image, (-1, self.noisy_image.shape[2]))
        full_cov = np.cov(reshaped, rowvar=False)

        cov = mnf.estimate_image_covariance(self.noisy_image)
        np.testing.assert_array_almost_equal(cov.cov(), full_cov, 3)

    def test_combine_two_identical_covariances_yields_original_covariance(self):
        cov_1 = mnf.estimate_image_covariance(self.noisy_image)
        cov_2 = mnf.estimate_image_covariance(self.noisy_image)

        cov_1.update_with_covariance(cov_2)
        np.testing.assert_array_almost_equal(cov_1.cov(), cov_2.cov())

    def test_orthonormalized_transform_yields_same_denoised_result_as_original_transform(self):
        denoised = mnf.denoise(self.forward, self.inverse, 4, self.noisy_image)

        forward_orth, inverse_orth, _ = mnf.mnf_transformation_matrix(self.image_covariance, self.noise_covariance, orthonormalize_bands=4)
        denoised_orth = mnf.denoise(forward_orth, inverse_orth, 4, self.noisy_image)

        np.testing.assert_array_almost_equal(denoised, denoised_orth)

    def test_keep_num_bands_argument_yields_the_first_n_bands_in_transform(self):
        forward, _ = mnf.transform_forward(self.forward, self.noisy_image)
        forward_5, _ = mnf.transform_forward(self.forward, self.noisy_image, 5)

        np.testing.assert_array_almost_equal(forward[:,:,0:5], forward_5)

    def test_mnf_class_yields_same_properties_as_individual_calls(self):
        mnf_class = mnf.mnf(self.noisy_image, 4, orthonormalize=False, verbose=False)
        np.testing.assert_array_almost_equal(mnf_class.cov.cov(), self.image_covariance.cov())
        np.testing.assert_array_almost_equal(mnf_class.noise_cov.cov(), self.noise_covariance.cov())
        np.testing.assert_array_almost_equal(mnf_class.forward, self.forward)
        np.testing.assert_array_almost_equal(mnf_class.inverse, self.inverse)

    def test_mnf_objects_from_parts_of_image_combines_to_mnf_transform_from_full_image(self):
        mnfs = []
        for i in np.arange(0, self.noisy_image.shape[0]):
            line = np.array([self.noisy_image[i, :, :]])
            transform = mnf.mnf(line, 4, False, False)
            mnfs.append(transform)

        for i in np.arange(1, len(mnfs)):
            mnfs[0].add_covariances(mnfs[i].cov, mnfs[i].noise_cov)

        np.testing.assert_array_almost_equal(mnfs[0].inverse, self.inverse)
        np.testing.assert_array_almost_equal(mnfs[0].forward, self.forward)
        np.testing.assert_array_almost_equal(mnfs[0].cov.cov(), self.image_covariance.cov())
        np.testing.assert_array_almost_equal(mnfs[0].noise_cov.cov(), self.noise_covariance.cov())

    def get_mnf_instance_and_two_spectra(self):
        """
        Get orthonormalized mnf instance and two denoised spectra and two
        spectra within MNF transform space for testing of dot product
        and euclidean distance.
        """
        mnf_class = mnf.mnf(self.noisy_image, 4, orthonormalize=True, verbose=False)

        c_1 = [57, 89]
        c_2 = [130, 140]

        denoised = mnf.denoise(self.forward, self.inverse, mnf_class.num_bands_in_inverse, self.noisy_image)
        spec_1 = denoised[c_1[0], c_1[1], :]
        spec_2 = denoised[c_2[0], c_2[1], :]

        transformed = mnf_class.transform(self.noisy_image)

        transformed_1 = transformed[c_1[0], c_1[1], :]
        transformed_2 = transformed[c_2[0], c_2[1], :]
        return spec_1, spec_2, transformed_1, transformed_2, mnf_class

    def test_euclidean_distances_in_orthonormalized_transform_are_the_same_as_distances_between_denoised_spectra(self):

        spec_1, spec_2, transformed_1, transformed_2, mnf_class = self.get_mnf_instance_and_two_spectra()

        self.assertAlmostEqual(np.sum((spec_1 - spec_2)**2), np.sum((transformed_1 - transformed_2)**2))

    def test_mnf_transformed_dot_product_equal_to_denoised_dot_product(self):
        spec_1, spec_2, transformed_1, transformed_2, mnf_class = self.get_mnf_instance_and_two_spectra()

        dot_product = np.sum(spec_1 * spec_2)

        transformed_dot_product = mnf_class.dot_product(transformed_1, transformed_2)

        self.assertAlmostEqual(dot_product, transformed_dot_product)

    def test_mnf_transformed_dot_product_for_multiple_spectra(self):
        _, _, transformed_1, transformed_2, mnf_class = self.get_mnf_instance_and_two_spectra()
        transformed_collection = np.array([transformed_2, transformed_2]).T
        transformed_dot_product = mnf_class.dot_product(transformed_1, transformed_2)
        transformed_dot_products = mnf_class.dot_product(transformed_1, transformed_collection)

        for dot in transformed_dot_products:
            self.assertAlmostEqual(dot, transformed_dot_product)

    def test_mnf_vector_length(self):
        spec_1, spec_2, transformed_1, transformed_2, mnf_class = self.get_mnf_instance_and_two_spectra()

        transformed_length = mnf_class.vector_length(np.array([[transformed_1]]))
        length = np.sqrt(np.sum(spec_1**2))

        self.assertAlmostEqual(transformed_length, length)




if __name__ == '__main__':
    unittest.main()
