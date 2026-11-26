# Copyright (c) d-SEAMS developers
"""ctypes calls for the pore geometry C API."""

from __future__ import annotations

import ctypes
from pathlib import Path

OK = 0
NULL = 1
RANGE = 2
BOX = 3
ARG = 4
STATE = 5


class PoreLib:
    def __init__(self, path: Path) -> None:
        self.lib = ctypes.CDLL(str(path))
        self.lib.seams_system_create.restype = ctypes.c_void_p
        self.lib.seams_system_destroy.argtypes = [ctypes.c_void_p]
        self.lib.seams_system_set_atoms.argtypes = [
            ctypes.c_void_p,
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_size_t,
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_size_t,
            ctypes.POINTER(ctypes.c_int),
            ctypes.c_void_p,
        ]
        self.lib.seams_system_set_atoms.restype = ctypes.c_int
        dptr = ctypes.POINTER(ctypes.c_double)
        iptr = ctypes.POINTER(ctypes.c_int)
        sptr = ctypes.POINTER(ctypes.c_size_t)
        self.lib.seams_density_fractional.argtypes = [
            ctypes.c_void_p,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_int,
            dptr,
            ctypes.c_size_t,
            ctypes.c_double,
            dptr,
            ctypes.c_size_t,
            dptr,
            dptr,
            dptr,
            dptr,
            ctypes.c_size_t,
            sptr,
        ]
        self.lib.seams_density_fractional.restype = ctypes.c_int
        self.lib.seams_pore_frame.argtypes = [
            ctypes.c_void_p,
            ctypes.c_int,
            ctypes.c_int,
            dptr,
            ctypes.c_size_t,
            dptr,
            ctypes.c_size_t,
            dptr,
            dptr,
            dptr,
            ctypes.POINTER(ctypes.c_int),
        ]
        self.lib.seams_pore_frame.restype = ctypes.c_int
        self.lib.seams_pore_profiles.argtypes = [
            ctypes.c_void_p,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_double,
            ctypes.c_int,
            dptr,
            ctypes.c_size_t,
            dptr,
            ctypes.c_size_t,
            ctypes.c_int,
            ctypes.c_double,
            dptr,
            ctypes.c_size_t,
            dptr,
            dptr,
            dptr,
            dptr,
            ctypes.c_size_t,
            sptr,
            dptr,
            dptr,
            dptr,
            dptr,
            ctypes.c_size_t,
            sptr,
            dptr,
            dptr,
            dptr,
            dptr,
            ctypes.c_size_t,
            sptr,
        ]
        self.lib.seams_pore_profiles.restype = ctypes.c_int
        self.lib.seams_host_guest.argtypes = [
            ctypes.c_void_p,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_double,
            iptr,
            iptr,
            iptr,
            dptr,
            ctypes.c_size_t,
            sptr,
            iptr,
            ctypes.c_size_t,
            sptr,
            iptr,
            ctypes.c_size_t,
        ]
        self.lib.seams_host_guest.restype = ctypes.c_int
        self.lib.seams_last_error.restype = ctypes.c_char_p

    def system(self, xyz, box, types):
        n = len(xyz) // 3
        handle = self.lib.seams_system_create()
        if not handle:
            raise RuntimeError("seams_system_create returned null")
        cxyz = (ctypes.c_double * len(xyz))(*xyz)
        cbox = (ctypes.c_double * 9)(*box)
        ctypes_types = (ctypes.c_int * n)(*types) if types is not None else None
        status = self.lib.seams_system_set_atoms(
            handle, cxyz, n, cbox, 9, ctypes_types, None
        )
        if status != OK:
            self.lib.seams_system_destroy(handle)
            raise RuntimeError(self.message())
        return handle

    def close(self, handle) -> None:
        self.lib.seams_system_destroy(handle)

    def message(self) -> str:
        text = self.lib.seams_last_error()
        return text.decode() if text else ""

    def fractional(self, handle, n, bins, axis, origin=0.0, direction=None, masses=None, type_id=0):
        coord = (ctypes.c_double * bins)()
        number = (ctypes.c_double * bins)()
        mass = (ctypes.c_double * bins)()
        volume = (ctypes.c_double * bins)()
        count = ctypes.c_size_t()
        direction_buf = None
        direction_len = 0
        if direction is not None:
            direction_buf = (ctypes.c_double * 3)(*direction)
            direction_len = 3
        mass_buf = None
        mass_len = 0
        if masses is not None:
            mass_buf = (ctypes.c_double * len(masses))(*masses)
            mass_len = len(masses)
        status = self.lib.seams_density_fractional(
            handle,
            type_id,
            bins,
            axis,
            direction_buf,
            direction_len,
            origin,
            mass_buf,
            mass_len,
            coord,
            number,
            mass,
            volume,
            bins,
            ctypes.byref(count),
        )
        if status != OK:
            raise RuntimeError(self.message())
        return {
            "coordinate": [coord[i] for i in range(count.value)],
            "number_density": [number[i] for i in range(count.value)],
            "mass_density": [mass[i] for i in range(count.value)],
            "volume": [volume[i] for i in range(count.value)],
        }

    def frame(self, handle, framework=0, axis=-1, direction=None, centre=None):
        axis_out = (ctypes.c_double * 3)()
        centre_out = (ctypes.c_double * 3)()
        frac_out = (ctypes.c_double * 3)()
        index = ctypes.c_int()
        direction_buf = None
        direction_len = 0
        if direction is not None:
            direction_buf = (ctypes.c_double * 3)(*direction)
            direction_len = 3
        centre_buf = None
        centre_len = 0
        if centre is not None:
            centre_buf = (ctypes.c_double * 3)(*centre)
            centre_len = 3
        status = self.lib.seams_pore_frame(
            handle,
            framework,
            axis,
            direction_buf,
            direction_len,
            centre_buf,
            centre_len,
            axis_out,
            centre_out,
            frac_out,
            ctypes.byref(index),
        )
        if status != OK:
            raise RuntimeError(self.message())
        return {
            "axis": [axis_out[i] for i in range(3)],
            "centre": [centre_out[i] for i in range(3)],
            "centre_frac": [frac_out[i] for i in range(3)],
            "axis_index": index.value,
        }

    def profiles(
        self,
        handle,
        radial_bins,
        azimuthal_bins,
        axial_bins,
        r_max,
        axis,
        centre,
        origin=0.0,
        type_id=0,
        framework=0,
    ):
        def column(n):
            return (ctypes.c_double * n)()

        radius = column(radial_bins)
        radial_number = column(radial_bins)
        radial_mass = column(radial_bins)
        radial_volume = column(radial_bins)
        theta = column(azimuthal_bins)
        azimuthal_number = column(azimuthal_bins)
        azimuthal_mass = column(azimuthal_bins)
        azimuthal_volume = column(azimuthal_bins)
        axial = column(axial_bins)
        axial_number = column(axial_bins)
        axial_mass = column(axial_bins)
        axial_volume = column(axial_bins)
        n_radial = ctypes.c_size_t()
        n_azimuthal = ctypes.c_size_t()
        n_axial = ctypes.c_size_t()
        centre_buf = (ctypes.c_double * 3)(*centre)
        status = self.lib.seams_pore_profiles(
            handle,
            type_id,
            radial_bins,
            azimuthal_bins,
            axial_bins,
            r_max,
            axis,
            None,
            0,
            centre_buf,
            3,
            framework,
            origin,
            None,
            0,
            radius,
            radial_number,
            radial_mass,
            radial_volume,
            radial_bins,
            ctypes.byref(n_radial),
            theta,
            azimuthal_number,
            azimuthal_mass,
            azimuthal_volume,
            azimuthal_bins,
            ctypes.byref(n_azimuthal),
            axial,
            axial_number,
            axial_mass,
            axial_volume,
            axial_bins,
            ctypes.byref(n_axial),
        )
        if status != OK:
            raise RuntimeError(self.message())

        def take(buf, n):
            return [buf[i] for i in range(n)]

        return {
            "radial": {
                "coordinate": take(radius, n_radial.value),
                "number_density": take(radial_number, n_radial.value),
                "mass_density": take(radial_mass, n_radial.value),
                "volume": take(radial_volume, n_radial.value),
            },
            "azimuthal": {
                "coordinate": take(theta, n_azimuthal.value),
                "number_density": take(azimuthal_number, n_azimuthal.value),
                "mass_density": take(azimuthal_mass, n_azimuthal.value),
                "volume": take(azimuthal_volume, n_azimuthal.value),
            },
            "axial": {
                "coordinate": take(axial, n_axial.value),
                "number_density": take(axial_number, n_axial.value),
                "mass_density": take(axial_mass, n_axial.value),
                "volume": take(axial_volume, n_axial.value),
            },
        }

    def host_guest(self, handle, guest, framework, cutoff):
        n_guest = ctypes.c_size_t()
        n_types = ctypes.c_size_t()
        status = self.lib.seams_host_guest(
            handle,
            guest,
            framework,
            cutoff,
            None,
            None,
            None,
            None,
            0,
            ctypes.byref(n_guest),
            None,
            0,
            ctypes.byref(n_types),
            None,
            0,
        )
        if status != OK:
            raise RuntimeError(self.message())
        guests = n_guest.value
        types = n_types.value
        guest_index = (ctypes.c_int * guests)()
        nearest = (ctypes.c_int * guests)()
        nearest_type = (ctypes.c_int * guests)()
        distance = (ctypes.c_double * max(guests, 1))()
        framework_types = (ctypes.c_int * max(types, 1))()
        contacts = (ctypes.c_int * max(guests * types, 1))()
        status = self.lib.seams_host_guest(
            handle,
            guest,
            framework,
            cutoff,
            guest_index,
            nearest,
            nearest_type,
            distance,
            guests,
            ctypes.byref(n_guest),
            framework_types,
            types,
            ctypes.byref(n_types),
            contacts,
            guests * types,
        )
        if status != OK:
            raise RuntimeError(self.message())
        return {
            "guest_index": [guest_index[i] for i in range(n_guest.value)],
            "nearest_index": [nearest[i] for i in range(n_guest.value)],
            "nearest_type": [nearest_type[i] for i in range(n_guest.value)],
            "distance": [distance[i] for i in range(n_guest.value)],
            "framework_types": [framework_types[i] for i in range(n_types.value)],
            "contacts": [
                contacts[g * n_types.value + t]
                for g in range(n_guest.value)
                for t in range(n_types.value)
            ],
        }

    def status_fractional(self, handle, bins, axis, direction=None, masses=None, n_atoms=0):
        count = ctypes.c_size_t()
        direction_buf = None
        direction_len = 0
        if direction is not None:
            direction_buf = (ctypes.c_double * len(direction))(*direction)
            direction_len = len(direction)
        mass_buf = None
        mass_len = 0
        if masses is not None:
            mass_buf = (ctypes.c_double * len(masses))(*masses)
            mass_len = len(masses)
        return self.lib.seams_density_fractional(
            handle,
            0,
            bins,
            axis,
            direction_buf,
            direction_len,
            0.0,
            mass_buf,
            mass_len,
            None,
            None,
            None,
            None,
            0,
            ctypes.byref(count),
        )
