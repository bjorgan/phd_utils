import unittest
import hyperspectral_utils
import tempfile
import numpy as np

class test_hyperread(unittest.TestCase):
    def setUp(self):
        #dummy data for writing and reading
        self.dummy_data = np.array([[[1, 2, 3, 4], [4, 5, 6, 7], [7, 8, 9, 10]],
                                    [[10, 11, 12, 13], [13, 14, 15, 16], [16, 17, 18, 20]]]).astype(np.float32)
        self.dummy_header = {'wlens': [67.0, 57.6, 89.0], 'default_bands': [0, 1, 2]}
        self.num_times_to_write = 7

    def test_writing_and_reading(self):
        temp_dir = tempfile.TemporaryDirectory()
        temp_dir_name = temp_dir.name

        #write image
        temp_image_base = temp_dir_name + '/image'
        hsi_file = hyperspectral_utils.hyperwrite(temp_image_base, self.dummy_header)
        for i in range(0, self.num_times_to_write):
            hsi_file.write(self.dummy_data)
        hsi_file.close()

        #read image
        hsi_image = hyperspectral_utils.hyperread(temp_image_base + '.img')

        #check that headers are equal
        self.assertEqual(hsi_image.header['lines'], self.dummy_data.shape[0]*self.num_times_to_write)
        self.assertEqual(hsi_image.header['samples'], self.dummy_data.shape[1])
        self.assertEqual(hsi_image.header['bands'], self.dummy_data.shape[2])
        for i, wlen in enumerate(hsi_image.header['wlens']):
            self.assertEqual(wlen, self.dummy_header['wlens'][i])
        for i, wlen in enumerate(hsi_image.header['default_bands']):
            self.assertEqual(wlen, self.dummy_header['default_bands'][i])

        #check that data are equal
        for i, data in enumerate(np.ravel(self.dummy_data)):
            self.assertEqual(data, np.ravel(hsi_image.image())[i])

        temp_dir.cleanup()

if __name__ == '__main__':
    unittest.main()
