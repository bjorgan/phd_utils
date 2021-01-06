Absorption spectra put here will be automatically loaded by
`forward_model.get_absorption_spectra()`. Spectra are assumed to be in
two-column format,

```
#wavelength  absorption
400  10
410  20
420  21
...
```

Use e.g. spectra from https://omlc.org/spectra/ or
https://omlc.org/spectra/PhotochemCAD/, and make sure the header is prepended
by `#`. (This repository includes only examples of blood absoption spectra.)
The filename is used as column header in the resulting absorption dataframe.
