# In-situ tutorial

Eight cubic-ice atoms, one observation through `InSitu.observe`.

```{eval-rst}
.. include:: ../videos/embeds/insitu.rst
```

```bash
meson setup build -Dwith_tests=true
meson compile -C build
python3 engines/insitu_observe.py build
```

`engines/insitu_observe.py` loads `libseams_lang.so` from the build tree and passes the eight-atom diamond with an orthorhombic box. It stops when NumPy or `libseams_lang.so` is absent. Otherwise every label is cubic and the script prints:

```
insitu atoms 8 cubic 8
```
