img(:,:,1) = reshape(1:20, [4 5])';
for i=2:20
	img(:,:,i) = 20*(i - 1) + img(:,:,1);
end
wavelengths = 20:40-1;

writeimage(single(img), wavelengths, 'data/test_image');
