hyperspectral_utils
===================

Python utilities for reading and writing hyperspectral images in ENVI format.

Installation instructions
-------------------------

To install locally, run

```
python3 setup.py install --user
```

numpy is required.

Usage example
-------------

Read image by

```
from hyperspectral_utils import hyperread
img_container = hyperread('path/to/image.hyspex')
```

Header information like wavelengths is accessed by e.g.
`img_container.header['wlens']`.  Contained image data is accessed by

```
image = img_container.image()
```

`image` is now a `NUM_LINES` (height) x `NUM_SAMPLES` (width) x
`NUM_WAVELENGTHS` shaped (memory mapped) numpy array.

See docstrings in `hyperspectral_utils.py` or `help(hyperread)` for details.
