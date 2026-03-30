#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to create v6-format HDF5 test files from existing v5 files.
These files use the new CADET v6 interface:
- UNIT_TYPE = COLUMN_MODEL_1D / COLUMN_MODEL_2D
- NPARTYPE at unit level (mandatory)
- Particle info in particle_type_xxx subgroups
- HAS_FILM_DIFFUSION, HAS_PORE_DIFFUSION, HAS_SURFACE_DIFFUSION flags
- PORE_DIFFUSION instead of PAR_DIFFUSION
- SURFACE_DIFFUSION instead of PAR_SURFDIFFUSION
- adsorption subgroup inside particle_type_xxx
"""

import h5py
import numpy as np
import os

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "data", "CADET_configs")


def create_base_structure(f, nunits=3):
    """Create the basic input/model structure."""
    f.create_dataset("input/model/NUNITS", data=nunits)


def create_inlet(f, unit_idx, ncomp=1):
    """Create an inlet unit."""
    grp = f.create_group(f"input/model/unit_{unit_idx:03d}")
    grp.create_dataset("UNIT_TYPE", data="INLET")
    grp.create_dataset("NCOMP", data=ncomp)
    return grp


def create_outlet(f, unit_idx, ncomp=1):
    """Create an outlet unit."""
    grp = f.create_group(f"input/model/unit_{unit_idx:03d}")
    grp.create_dataset("UNIT_TYPE", data="OUTLET")
    grp.create_dataset("NCOMP", data=ncomp)
    return grp


def create_v6_PlugFlow():
    """Plug flow: COLUMN_MODEL_1D, NPARTYPE=0, no particles."""
    fname = os.path.join(OUTPUT_DIR, "v6_PlugFlow_1comp.h5")
    with h5py.File(fname, 'w') as f:
        create_base_structure(f)
        create_inlet(f, 0)

        unit = f.create_group("input/model/unit_001")
        unit.create_dataset("UNIT_TYPE", data="COLUMN_MODEL_1D")
        unit.create_dataset("NCOMP", data=1)
        unit.create_dataset("NPARTYPE", data=0)
        unit.create_dataset("COL_DISPERSION", data=0.0)
        unit.create_dataset("COL_LENGTH", data=1.0)
        unit.create_dataset("TOTAL_POROSITY", data=1.0)
        unit.create_dataset("VELOCITY", data=0.03333333333333333)
        unit.create_dataset("INIT_C", data=np.array([0.0]))

        disc = unit.create_group("discretization")
        disc.create_dataset("NELEM", data=1)
        disc.create_dataset("POLYDEG", data=3)
        disc.create_dataset("SPATIAL_METHOD", data="DG")

        create_outlet(f, 2)
    print(f"Created {fname}")


def create_v6_LRM():
    """LRM: COLUMN_MODEL_1D with equilibrium particle (no film/pore/surface diffusion)."""
    fname = os.path.join(OUTPUT_DIR, "v6_LRM_dynLin_1comp.h5")
    with h5py.File(fname, 'w') as f:
        create_base_structure(f)
        create_inlet(f, 0)

        unit = f.create_group("input/model/unit_001")
        unit.create_dataset("UNIT_TYPE", data="COLUMN_MODEL_1D")
        unit.create_dataset("NCOMP", data=1)
        unit.create_dataset("NPARTYPE", data=1)
        unit.create_dataset("COL_DISPERSION", data=0.0001)
        unit.create_dataset("COL_LENGTH", data=1.0)
        unit.create_dataset("TOTAL_POROSITY", data=0.6)
        unit.create_dataset("VELOCITY", data=0.03333333333333333)
        unit.create_dataset("INIT_C", data=np.array([0.0]))
        unit.create_dataset("INIT_Q", data=np.array([0.0]))

        disc = unit.create_group("discretization")
        disc.create_dataset("NELEM", data=1)
        disc.create_dataset("POLYDEG", data=3)
        disc.create_dataset("SPATIAL_METHOD", data="DG")

        # particle_type_000: equilibrium particle (LRM-like)
        pt = unit.create_group("particle_type_000")
        pt.create_dataset("HAS_FILM_DIFFUSION", data=False)
        pt.create_dataset("HAS_PORE_DIFFUSION", data=False)
        pt.create_dataset("HAS_SURFACE_DIFFUSION", data=False)
        pt.create_dataset("ADSORPTION_MODEL", data="LINEAR")
        pt.create_dataset("NBOUND", data=np.array([1]))

        ads = pt.create_group("adsorption")
        ads.create_dataset("IS_KINETIC", data=1)
        ads.create_dataset("LIN_KA", data=np.array([1.0]))
        ads.create_dataset("LIN_KD", data=np.array([1.0]))

        create_outlet(f, 2)
    print(f"Created {fname}")


def create_v6_LRMP():
    """LRMP: COLUMN_MODEL_1D with film diffusion only (homogeneous particle)."""
    fname = os.path.join(OUTPUT_DIR, "v6_LRMP_dynLin_1comp.h5")
    with h5py.File(fname, 'w') as f:
        create_base_structure(f)
        create_inlet(f, 0)

        unit = f.create_group("input/model/unit_001")
        unit.create_dataset("UNIT_TYPE", data="COLUMN_MODEL_1D")
        unit.create_dataset("NCOMP", data=1)
        unit.create_dataset("NPARTYPE", data=1)
        unit.create_dataset("COL_DISPERSION", data=np.array([0.0001]))
        unit.create_dataset("COL_LENGTH", data=1.0)
        unit.create_dataset("COL_POROSITY", data=0.6)
        unit.create_dataset("VELOCITY", data=0.03333333333333333)
        unit.create_dataset("INIT_C", data=np.array([0.0]))
        unit.create_dataset("INIT_CP", data=np.array([0.0]))
        unit.create_dataset("INIT_Q", data=np.array([0.0]))
        unit.create_dataset("PAR_TYPE_VOLFRAC", data=1)

        disc = unit.create_group("discretization")
        disc.create_dataset("NELEM", data=1)
        disc.create_dataset("POLYDEG", data=3)
        disc.create_dataset("SPATIAL_METHOD", data="DG")
        disc.create_dataset("PAR_GEOM", data=np.array([b'SPHERE']))

        # particle_type_000: homogeneous particle (LRMP-like)
        pt = unit.create_group("particle_type_000")
        pt.create_dataset("HAS_FILM_DIFFUSION", data=True)
        pt.create_dataset("HAS_PORE_DIFFUSION", data=False)
        pt.create_dataset("HAS_SURFACE_DIFFUSION", data=False)
        pt.create_dataset("PAR_POROSITY", data=0.2)
        pt.create_dataset("PAR_RADIUS", data=0.0001)
        pt.create_dataset("FILM_DIFFUSION", data=np.array([0.00333333]))
        pt.create_dataset("ADSORPTION_MODEL", data="LINEAR")
        pt.create_dataset("NBOUND", data=np.array([1]))

        ads = pt.create_group("adsorption")
        ads.create_dataset("IS_KINETIC", data=1)
        ads.create_dataset("LIN_KA", data=np.array([1.0]))
        ads.create_dataset("LIN_KD", data=np.array([1.0]))

        create_outlet(f, 2)
    print(f"Created {fname}")


def create_v6_GRM():
    """GRM: COLUMN_MODEL_1D with film + pore diffusion, no surface diffusion."""
    fname = os.path.join(OUTPUT_DIR, "v6_GRM_dynLin_1comp.h5")
    with h5py.File(fname, 'w') as f:
        create_base_structure(f)
        create_inlet(f, 0)

        unit = f.create_group("input/model/unit_001")
        unit.create_dataset("UNIT_TYPE", data="COLUMN_MODEL_1D")
        unit.create_dataset("NCOMP", data=1)
        unit.create_dataset("NPARTYPE", data=1)
        unit.create_dataset("COL_DISPERSION", data=5.75e-08)
        unit.create_dataset("COL_LENGTH", data=0.014)
        unit.create_dataset("COL_POROSITY", data=0.37)
        unit.create_dataset("VELOCITY", data=0.000575)
        unit.create_dataset("INIT_C", data=np.array([0.0]))
        unit.create_dataset("INIT_CP", data=np.array([0.0]))
        unit.create_dataset("INIT_Q", data=np.array([0.0]))
        unit.create_dataset("PAR_TYPE_VOLFRAC", data=1)

        disc = unit.create_group("discretization")
        disc.create_dataset("NELEM", data=8)
        disc.create_dataset("POLYDEG", data=3)
        disc.create_dataset("SPATIAL_METHOD", data="DG")
        disc.create_dataset("PAR_GEOM", data=np.array([b'SPHERE']))

        # particle_type_000: GRM particle (film + pore diffusion)
        pt = unit.create_group("particle_type_000")
        pt.create_dataset("HAS_FILM_DIFFUSION", data=True)
        pt.create_dataset("HAS_PORE_DIFFUSION", data=True)
        pt.create_dataset("HAS_SURFACE_DIFFUSION", data=False)
        pt.create_dataset("PAR_POROSITY", data=0.75)
        pt.create_dataset("PAR_RADIUS", data=4.5e-05)
        pt.create_dataset("PAR_CORERADIUS", data=0.0)
        pt.create_dataset("FILM_DIFFUSION", data=np.array([6.9e-06]))
        pt.create_dataset("PORE_DIFFUSION", data=np.array([6.07e-11]))
        pt.create_dataset("ADSORPTION_MODEL", data="LINEAR")
        pt.create_dataset("NBOUND", data=np.array([1]))

        ads = pt.create_group("adsorption")
        ads.create_dataset("IS_KINETIC", data=True)
        ads.create_dataset("LIN_KA", data=np.array([3.55]))
        ads.create_dataset("LIN_KD", data=np.array([0.1]))

        create_outlet(f, 2)
    print(f"Created {fname}")


def create_v6_GRMsd():
    """GRMsd: COLUMN_MODEL_1D with film + pore + surface diffusion."""
    fname = os.path.join(OUTPUT_DIR, "v6_GRMsd_dynLin_1comp.h5")
    with h5py.File(fname, 'w') as f:
        create_base_structure(f)
        create_inlet(f, 0)

        unit = f.create_group("input/model/unit_001")
        unit.create_dataset("UNIT_TYPE", data="COLUMN_MODEL_1D")
        unit.create_dataset("NCOMP", data=1)
        unit.create_dataset("NPARTYPE", data=1)
        unit.create_dataset("COL_DISPERSION", data=5.75e-08)
        unit.create_dataset("COL_LENGTH", data=0.014)
        unit.create_dataset("COL_POROSITY", data=0.37)
        unit.create_dataset("VELOCITY", data=0.000575)
        unit.create_dataset("INIT_C", data=np.array([0.0]))
        unit.create_dataset("INIT_CP", data=np.array([0.0]))
        unit.create_dataset("INIT_Q", data=np.array([0.0]))
        unit.create_dataset("PAR_TYPE_VOLFRAC", data=1)

        disc = unit.create_group("discretization")
        disc.create_dataset("NELEM", data=8)
        disc.create_dataset("POLYDEG", data=3)
        disc.create_dataset("SPATIAL_METHOD", data="DG")
        disc.create_dataset("PAR_GEOM", data=np.array([b'SPHERE']))

        # particle_type_000: GRM particle with surface diffusion
        pt = unit.create_group("particle_type_000")
        pt.create_dataset("HAS_FILM_DIFFUSION", data=True)
        pt.create_dataset("HAS_PORE_DIFFUSION", data=True)
        pt.create_dataset("HAS_SURFACE_DIFFUSION", data=True)
        pt.create_dataset("PAR_POROSITY", data=0.75)
        pt.create_dataset("PAR_RADIUS", data=4.5e-05)
        pt.create_dataset("PAR_CORERADIUS", data=0.0)
        pt.create_dataset("FILM_DIFFUSION", data=np.array([6.9e-06]))
        pt.create_dataset("PORE_DIFFUSION", data=np.array([6.07e-11]))
        pt.create_dataset("SURFACE_DIFFUSION", data=np.array([6.07e-11]))
        pt.create_dataset("ADSORPTION_MODEL", data="LINEAR")
        pt.create_dataset("NBOUND", data=np.array([1]))

        ads = pt.create_group("adsorption")
        ads.create_dataset("IS_KINETIC", data=True)
        ads.create_dataset("LIN_KA", data=np.array([3.55]))
        ads.create_dataset("LIN_KD", data=np.array([0.1]))

        create_outlet(f, 2)
    print(f"Created {fname}")


def create_v6_GRMsd_PSD():
    """GRMsd with 2 particle types (PSD)."""
    fname = os.path.join(OUTPUT_DIR, "v6_GRMsd_PSD_dynLin_1comp.h5")
    with h5py.File(fname, 'w') as f:
        create_base_structure(f)
        create_inlet(f, 0)

        unit = f.create_group("input/model/unit_001")
        unit.create_dataset("UNIT_TYPE", data="COLUMN_MODEL_1D")
        unit.create_dataset("NCOMP", data=1)
        unit.create_dataset("NPARTYPE", data=2)
        unit.create_dataset("COL_DISPERSION", data=5.75e-08)
        unit.create_dataset("COL_LENGTH", data=0.014)
        unit.create_dataset("COL_POROSITY", data=0.37)
        unit.create_dataset("VELOCITY", data=0.000575)
        unit.create_dataset("INIT_C", data=np.array([0.0]))
        unit.create_dataset("INIT_CP", data=np.array([0.0]))
        unit.create_dataset("INIT_Q", data=np.array([0.0]))
        unit.create_dataset("PAR_TYPE_VOLFRAC", data=1)

        disc = unit.create_group("discretization")
        disc.create_dataset("NELEM", data=8)
        disc.create_dataset("POLYDEG", data=3)
        disc.create_dataset("SPATIAL_METHOD", data="DG")
        disc.create_dataset("PAR_GEOM", data=np.array([b'SPHERE']))

        # particle_type_000
        pt0 = unit.create_group("particle_type_000")
        pt0.create_dataset("HAS_FILM_DIFFUSION", data=True)
        pt0.create_dataset("HAS_PORE_DIFFUSION", data=True)
        pt0.create_dataset("HAS_SURFACE_DIFFUSION", data=True)
        pt0.create_dataset("PAR_POROSITY", data=0.75)
        pt0.create_dataset("PAR_RADIUS", data=4.5e-05)
        pt0.create_dataset("PAR_CORERADIUS", data=0.0)
        pt0.create_dataset("FILM_DIFFUSION", data=np.array([6.9e-06]))
        pt0.create_dataset("PORE_DIFFUSION", data=np.array([6.07e-11]))
        pt0.create_dataset("SURFACE_DIFFUSION", data=np.array([1e-11]))
        pt0.create_dataset("ADSORPTION_MODEL", data="LINEAR")
        pt0.create_dataset("NBOUND", data=np.array([1]))

        ads0 = pt0.create_group("adsorption")
        ads0.create_dataset("IS_KINETIC", data=True)
        ads0.create_dataset("LIN_KA", data=np.array([3.55]))
        ads0.create_dataset("LIN_KD", data=np.array([0.1]))

        # particle_type_001
        pt1 = unit.create_group("particle_type_001")
        pt1.create_dataset("HAS_FILM_DIFFUSION", data=True)
        pt1.create_dataset("HAS_PORE_DIFFUSION", data=True)
        pt1.create_dataset("HAS_SURFACE_DIFFUSION", data=True)
        pt1.create_dataset("PAR_POROSITY", data=0.75)
        pt1.create_dataset("PAR_RADIUS", data=4.5e-05)
        pt1.create_dataset("PAR_CORERADIUS", data=0.0)
        pt1.create_dataset("FILM_DIFFUSION", data=np.array([6.9e-06]))
        pt1.create_dataset("PORE_DIFFUSION", data=np.array([6.07e-11]))
        pt1.create_dataset("SURFACE_DIFFUSION", data=np.array([1e-11]))
        pt1.create_dataset("ADSORPTION_MODEL", data="LINEAR")
        pt1.create_dataset("NBOUND", data=np.array([1]))

        ads1 = pt1.create_group("adsorption")
        ads1.create_dataset("IS_KINETIC", data=True)
        ads1.create_dataset("LIN_KA", data=np.array([5.0]))
        ads1.create_dataset("LIN_KD", data=np.array([0.1]))

        create_outlet(f, 2)
    print(f"Created {fname}")


def create_v6_GRMsd2D():
    """GRMsd 2D: COLUMN_MODEL_2D with film + pore + surface diffusion.
    Mirrors the v5 2DGRMsd3Zone file: unit_000=column, units 1-3=inlets, units 4-6=outlets."""
    fname = os.path.join(OUTPUT_DIR, "v6_GRMsd2D_dynLin_1comp.h5")
    with h5py.File(fname, 'w') as f:
        f.create_dataset("input/model/NUNITS", data=7)

        unit = f.create_group("input/model/unit_000")
        unit.create_dataset("UNIT_TYPE", data="COLUMN_MODEL_2D")
        unit.create_dataset("NCOMP", data=1)
        unit.create_dataset("NPARTYPE", data=1)
        unit.create_dataset("COL_DISPERSION", data=5.75e-08)
        unit.create_dataset("COL_DISPERSION_RADIAL", data=5e-08)
        unit.create_dataset("COL_LENGTH", data=0.014)
        unit.create_dataset("COL_POROSITY", data=0.37)
        unit.create_dataset("COL_RADIUS", data=0.0035)
        unit.create_dataset("CROSS_SECTION_AREA", data=3.848451000647497e-05)
        unit.create_dataset("VELOCITY", data=0.000575)
        unit.create_dataset("PORTS", data=3)
        unit.create_dataset("INIT_C", data=np.array([0]))
        unit.create_dataset("INIT_CP", data=np.array([0]))
        unit.create_dataset("INIT_Q", data=np.array([0]))
        unit.create_dataset("PAR_TYPE_VOLFRAC", data=1.0)

        disc = unit.create_group("discretization")
        disc.create_dataset("NCOL", data=4)
        disc.create_dataset("NPAR", data=3)
        disc.create_dataset("NRAD", data=3)
        disc.create_dataset("SPATIAL_METHOD", data="FV")
        disc.create_dataset("PAR_DISC_TYPE", data=np.array([b'EQUIDISTANT_PAR']))
        disc.create_dataset("RADIAL_DISC_TYPE", data="EQUIDISTANT")

        # particle_type_000: GRM with surface diffusion
        pt = unit.create_group("particle_type_000")
        pt.create_dataset("HAS_FILM_DIFFUSION", data=True)
        pt.create_dataset("HAS_PORE_DIFFUSION", data=True)
        pt.create_dataset("HAS_SURFACE_DIFFUSION", data=True)
        pt.create_dataset("PAR_POROSITY", data=0.75)
        pt.create_dataset("PAR_RADIUS", data=4.5e-05)
        pt.create_dataset("FILM_DIFFUSION", data=np.array([6.9e-06]))
        pt.create_dataset("PORE_DIFFUSION", data=np.array([6.07e-11]))
        pt.create_dataset("SURFACE_DIFFUSION", data=np.array([1e-11]))
        pt.create_dataset("NBOUND", data=1)
        pt.create_dataset("ADSORPTION_MODEL", data="LINEAR")

        ads = pt.create_group("adsorption")
        ads.create_dataset("IS_KINETIC", data=1)
        ads.create_dataset("LIN_KA", data=35.5)
        ads.create_dataset("LIN_KD", data=1.0)

        # Inlets: units 1-3
        for i in range(1, 4):
            grp = f.create_group(f"input/model/unit_{i:03d}")
            grp.create_dataset("UNIT_TYPE", data="INLET")
            grp.create_dataset("NCOMP", data=1)

        # Outlets: units 4-6
        for i in range(4, 7):
            grp = f.create_group(f"input/model/unit_{i:03d}")
            grp.create_dataset("UNIT_TYPE", data="OUTLET")
            grp.create_dataset("NCOMP", data=1)

    print(f"Created {fname}")


def create_v6_CSTR():
    """CSTR: stays as CSTR unit type, but with NPARTYPE=0 and no particles.
    Note: CSTR unit type name doesn't change in v6."""
    fname = os.path.join(OUTPUT_DIR, "v6_CSTR.h5")
    with h5py.File(fname, 'w') as f:
        create_base_structure(f)
        create_inlet(f, 0, ncomp=384)

        unit = f.create_group("input/model/unit_001")
        unit.create_dataset("UNIT_TYPE", data="CSTR")
        unit.create_dataset("NCOMP", data=384)
        unit.create_dataset("CONST_SOLID_VOLUME", data=0.0)
        unit.create_dataset("INIT_VOLUME", data=0.0005)
        unit.create_dataset("ADSORPTION_MODEL", data="NONE")
        unit.create_dataset("REACTION_MODEL", data="CRYSTALLIZATION")
        unit.create_dataset("USE_ANALYTIC_JACOBIAN", data=1)
        # Use a simple init array
        unit.create_dataset("INIT_C", data=np.zeros(384))

        create_outlet(f, 2, ncomp=384)
    print(f"Created {fname}")


if __name__ == "__main__":
    create_v6_PlugFlow()
    create_v6_LRM()
    create_v6_LRMP()
    create_v6_GRM()
    create_v6_GRMsd()
    create_v6_GRMsd_PSD()
    create_v6_GRMsd2D()
    create_v6_CSTR()
    print("\nAll v6 HDF5 test files created successfully!")
