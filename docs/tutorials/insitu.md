# In-situ tutorial

d-SEAMS 2.14.0. Eight cubic-ice atoms, one observation:

```bash
meson setup build -Dwith_tests=true
meson compile -C build
python3 engines/insitu_observe.py build
```
