import numpy as np

num_bands = 20
num_lines = 5
num_samples = 4

#generate image reflecting assumed input in readimage-t
img = np.zeros((num_lines, num_samples, num_bands), dtype=np.float32)
img[:,:,0] = np.arange(1, num_bands+1).reshape((num_lines, num_samples))
for i in range(1, num_bands):
	img[:,:,i] = 20*i + img[:,:,0]

#write image to file with BIL interleave
fid = open("test_image.img", "wb")
for i in range(num_lines):
    line_bil = img[i,:,:].swapaxes(-2, -1)
    line_bil.tofile(fid)
fid.close()

