# -*- coding: utf-8 -*-
"""
Unit tests for src/load_CADET_h5.py.
Each public function is tested in isolation using mock HDF5 groups.
"""

import pytest
import numpy as np
from unittest.mock import MagicMock
from src.load_CADET_h5 import (
    get_h5_value,
    map_unit_type_to_column_geometry,
    map_unit_type_to_column_model,
    map_unit_to_particle_model,
    is_v6_interface,
    extract_config_data_from_unit,
    extract_crystallization_config,
    _is_crystallization_unit,
    _get_crystallization_reaction_group,
    _extract_reaction_config,
    CADET_column_unit_types,
)


# %% Helpers

def _make_h5_group(mapping: dict, subgroups: dict = None):
    """Create a mock HDF5 group that mimics h5py dataset access.

    Keys present in *mapping* return mock datasets whose [()] call
    yields the corresponding value.  Missing keys return None.

    *subgroups* is an optional dict of {name: _make_h5_group(...)},
    accessible via group['name'] and the ``in`` operator.
    """
    subgroups = subgroups or {}
    group = MagicMock()

    def get_side_effect(k):
        if k in mapping:
            ds = MagicMock()
            val = mapping[k]
            type(ds).__getitem__ = lambda self, idx: val
            return ds
        if k in subgroups:
            return subgroups[k]
        return None

    group.get = get_side_effect
    group.__contains__ = lambda self, k: k in mapping or k in subgroups
    group.__getitem__ = lambda self, k: subgroups[k] if k in subgroups else get_side_effect(k)
    return group


# %% get_h5_value

@pytest.mark.ci
@pytest.mark.unit_test
def test_get_h5_value_missing_key():
    """A missing key should return None."""
    group = MagicMock()
    group.get.return_value = None
    assert get_h5_value(group, 'MISSING_KEY') is None


@pytest.mark.ci
@pytest.mark.unit_test
def test_get_h5_value_decodes_bytes():
    """Byte-string values should be decoded to str."""
    group = _make_h5_group({'UNIT_TYPE': b'GENERAL_RATE_MODEL'})
    assert get_h5_value(group, 'UNIT_TYPE') == 'GENERAL_RATE_MODEL'


@pytest.mark.ci
@pytest.mark.unit_test
def test_get_h5_value_scalar():
    """Scalar numeric values should be returned as-is."""
    group = _make_h5_group({'COL_POROSITY': 0.37})
    assert get_h5_value(group, 'COL_POROSITY') == 0.37


@pytest.mark.ci
@pytest.mark.unit_test
def test_get_h5_value_array_first_entry():
    """With firstEntryIfList=True (default), only the first element should be returned."""
    group = _make_h5_group({'COL_DISPERSION': np.array([5.75e-8, 1.0e-7])})
    assert get_h5_value(group, 'COL_DISPERSION', firstEntryIfList=True) == 5.75e-8


@pytest.mark.ci
@pytest.mark.unit_test
def test_get_h5_value_array_full():
    """With firstEntryIfList=False, the full array should be returned."""
    arr = np.array([5.75e-8, 1.0e-7])
    group = _make_h5_group({'COL_DISPERSION': arr})
    result = get_h5_value(group, 'COL_DISPERSION', firstEntryIfList=False)
    assert hasattr(result, '__len__')
    assert len(result) == 2


@pytest.mark.ci
@pytest.mark.unit_test
def test_get_h5_value_empty_array():
    """An empty array with firstEntryIfList=True should return the empty array (no crash)."""
    group = _make_h5_group({'EMPTY': np.array([])})
    result = get_h5_value(group, 'EMPTY', firstEntryIfList=True)
    assert hasattr(result, '__len__')
    assert len(result) == 0


# %% map_unit_type_to_column_model

@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("unit_type, expected", [
    ('GENERAL_RATE_MODEL',                        "Axial flow cylinder"),
    ('LUMPED_RATE_MODEL_WITHOUT_PORES',            "Axial flow cylinder"),
    ('CSTR',                                       "Axial flow cylinder"),
    ('COLUMN_MODEL_1D',                            "Axial flow cylinder"),
    ('RADIAL_GENERAL_RATE_MODEL',                  "Radial flow cylinder"),
    ('RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES',     "Radial flow cylinder"),
    ('RADIAL_COLUMN_MODEL_1D',                     "Radial flow cylinder"),
    ('FRUSTUM_GENERAL_RATE_MODEL',                 "Frustum"),
    ('FRUSTUM_LUMPED_RATE_MODEL_WITHOUT_PORES',    "Frustum"),
    ('FRUSTUM_COLUMN_MODEL_1D',                    "Frustum"),
])
def test_map_unit_type_to_column_geometry(unit_type, expected):
    """Every supported CADET unit type should map to its correct column geometry."""
    assert map_unit_type_to_column_geometry(unit_type) == expected


@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("unit_type, expected", [
    ('GENERAL_RATE_MODEL',                 "1D (axial coordinate)"),
    ('GENERAL_RATE_MODEL_DG',              "1D (axial coordinate)"),
    ('LUMPED_RATE_MODEL_WITHOUT_PORES',    "1D (axial coordinate)"),
    ('LUMPED_RATE_MODEL_WITH_PORES',       "1D (axial coordinate)"),
    ('LUMPED_RATE_MODEL_WITHOUT_PORES_DG', "1D (axial coordinate)"),
    ('LUMPED_RATE_MODEL_WITH_PORES_DG',    "1D (axial coordinate)"),
    ('FRUSTUM_GENERAL_RATE_MODEL',         "1D (axial coordinate)"),
    ('FRUSTUM_COLUMN_MODEL_1D',            "1D (axial coordinate)"),
    ('GENERAL_RATE_MODEL_2D',              "2D (axial and radial coordinate)"),
    ('CSTR',                               "0D (Homogeneous Tank)"),
    ('RADIAL_GENERAL_RATE_MODEL',          "1D (radial coordinate)"),
    ('RADIAL_COLUMN_MODEL_1D',             "1D (radial coordinate)"),
])
def test_map_unit_type_to_column_model(unit_type, expected):
    """Every supported CADET unit type should map to its correct column resolution."""
    assert map_unit_type_to_column_model(unit_type) == expected


@pytest.mark.ci
@pytest.mark.unit_test
def test_map_unit_type_to_column_model_invalid():
    """An unsupported unit type should raise ValueError."""
    with pytest.raises(ValueError, match="Invalid unit type"):
        map_unit_type_to_column_model('INVALID_MODEL')


# %% map_unit_to_particle_model

@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("unit_type, h5_values, expected", [
    ('GENERAL_RATE_MODEL',              {'COL_POROSITY': 0.37},                                 "1D (radial coordinate)"),
    ('GENERAL_RATE_MODEL',              {'COL_POROSITY': 1.0},                                  None),
    ('LUMPED_RATE_MODEL_WITH_PORES',    {'COL_POROSITY': 0.37},                                 "0D (homogeneous)"),
    ('LUMPED_RATE_MODEL_WITHOUT_PORES', {'TOTAL_POROSITY': 1.0},                                None),
    ('LUMPED_RATE_MODEL_WITHOUT_PORES', {'TOTAL_POROSITY': 0.5, 'CONST_SOLID_VOLUME': 0.0},     None),
    ('CSTR',                            {'TOTAL_POROSITY': 0.5, 'CONST_SOLID_VOLUME': 0.1},     "0D (homogeneous)"),
    ('CSTR',                            {'TOTAL_POROSITY': 1.0},                                None),
])
def test_map_unit_to_particle_model(unit_type, h5_values, expected):
    """Particle model mapping should depend on unit type and porosity/solid-volume values."""
    group = _make_h5_group(h5_values)
    assert map_unit_to_particle_model(unit_type, group) == expected


@pytest.mark.ci
@pytest.mark.unit_test
def test_map_unit_to_particle_model_invalid():
    """An unsupported unit type should raise ValueError."""
    group = _make_h5_group({})
    with pytest.raises(ValueError, match="Invalid unit type"):
        map_unit_to_particle_model('INVALID_MODEL', group)


# %% CADET_column_unit_types constant

@pytest.mark.ci
@pytest.mark.unit_test
def test_CADET_column_unit_types_completeness():
    """The constant should contain all supported CADET column unit types (v5 + v6)."""
    expected = [
        'CSTR',
        'COLUMN_MODEL_1D', 'COLUMN_MODEL_2D',
        'GENERAL_RATE_MODEL', 'LUMPED_RATE_MODEL_WITHOUT_PORES',
        'LUMPED_RATE_MODEL_WITH_PORES', 'GENERAL_RATE_MODEL_2D',
        'FRUSTUM_COLUMN_MODEL_1D',
        'FRUSTUM_GENERAL_RATE_MODEL', 'FRUSTUM_LUMPED_RATE_MODEL_WITHOUT_PORES',
        'FRUSTUM_LUMPED_RATE_MODEL_WITH_PORES',
        'RADIAL_COLUMN_MODEL_1D',
        'RADIAL_GENERAL_RATE_MODEL', 'RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES',
        'RADIAL_LUMPED_RATE_MODEL_WITH_PORES',
        'GENERAL_RATE_MODEL_DG', 'LUMPED_RATE_MODEL_WITHOUT_PORES_DG',
        'LUMPED_RATE_MODEL_WITH_PORES_DG',
    ]
    assert len(CADET_column_unit_types) == len(expected)
    for t in expected:
        assert t in CADET_column_unit_types


# %% is_v6_interface

@pytest.mark.ci
@pytest.mark.unit_test
def test_is_v6_interface_column_model_1d():
    """COLUMN_MODEL_1D should always be detected as v6."""
    group = _make_h5_group({})
    assert is_v6_interface('COLUMN_MODEL_1D', group) is True


@pytest.mark.ci
@pytest.mark.unit_test
def test_is_v6_interface_column_model_2d():
    """COLUMN_MODEL_2D should always be detected as v6."""
    group = _make_h5_group({})
    assert is_v6_interface('COLUMN_MODEL_2D', group) is True


@pytest.mark.ci
@pytest.mark.unit_test
def test_is_v6_interface_old_type_with_particle_type_group():
    """Old unit type with particle_type_000 subgroup should be detected as v6."""
    pt_group = _make_h5_group({'HAS_FILM_DIFFUSION': True})
    group = _make_h5_group({'NPARTYPE': 1}, subgroups={'particle_type_000': pt_group})
    assert is_v6_interface('GENERAL_RATE_MODEL', group) is True


@pytest.mark.ci
@pytest.mark.unit_test
def test_is_v6_interface_old_type_without_particle_type_group():
    """Old unit type without particle_type_000 is v5."""
    group = _make_h5_group({'NPARTYPE': 1})
    assert is_v6_interface('GENERAL_RATE_MODEL', group) is False


@pytest.mark.ci
@pytest.mark.unit_test
def test_is_v6_interface_old_type_npartype_zero():
    """Old unit type with NPARTYPE=0 is v5 (no interface changes for npartype=0)."""
    group = _make_h5_group({'NPARTYPE': 0})
    assert is_v6_interface('LUMPED_RATE_MODEL_WITHOUT_PORES', group) is False


# %% map_unit_to_particle_model – v6 paths

@pytest.mark.ci
@pytest.mark.unit_test
def test_map_v6_particle_model_column_model_1d_with_pore_diff():
    """COLUMN_MODEL_1D with HAS_PORE_DIFFUSION should map to 1D (radial coordinate)."""
    pt_group = _make_h5_group({
        'HAS_FILM_DIFFUSION': True,
        'HAS_PORE_DIFFUSION': True,
        'HAS_SURFACE_DIFFUSION': False,
    })
    group = _make_h5_group({'NPARTYPE': 1}, subgroups={'particle_type_000': pt_group})
    assert map_unit_to_particle_model('COLUMN_MODEL_1D', group) == "1D (radial coordinate)"


@pytest.mark.ci
@pytest.mark.unit_test
def test_map_v6_particle_model_film_diff_only():
    """COLUMN_MODEL_1D with only HAS_FILM_DIFFUSION should map to 0D (homogeneous)."""
    pt_group = _make_h5_group({
        'HAS_FILM_DIFFUSION': True,
        'HAS_PORE_DIFFUSION': False,
        'HAS_SURFACE_DIFFUSION': False,
    })
    group = _make_h5_group({'NPARTYPE': 1}, subgroups={'particle_type_000': pt_group})
    assert map_unit_to_particle_model('COLUMN_MODEL_1D', group) == "0D (homogeneous)"


@pytest.mark.ci
@pytest.mark.unit_test
def test_map_v6_particle_model_equilibrium():
    """COLUMN_MODEL_1D with no diffusion flags should map to 0D (homogeneous)."""
    pt_group = _make_h5_group({
        'HAS_FILM_DIFFUSION': False,
        'HAS_PORE_DIFFUSION': False,
        'HAS_SURFACE_DIFFUSION': False,
    })
    group = _make_h5_group({'NPARTYPE': 1}, subgroups={'particle_type_000': pt_group})
    assert map_unit_to_particle_model('COLUMN_MODEL_1D', group) == "0D (homogeneous)"


@pytest.mark.ci
@pytest.mark.unit_test
def test_map_v6_particle_model_npartype_zero():
    """COLUMN_MODEL_1D with NPARTYPE=0 should return None (no particles)."""
    group = _make_h5_group({'NPARTYPE': 0})
    assert map_unit_to_particle_model('COLUMN_MODEL_1D', group) is None


@pytest.mark.ci
@pytest.mark.unit_test
def test_map_v6_particle_model_missing_particle_type_group():
    """COLUMN_MODEL_1D with NPARTYPE=1 but no particle_type_000 should return None."""
    group = _make_h5_group({'NPARTYPE': 1})
    assert map_unit_to_particle_model('COLUMN_MODEL_1D', group) is None


@pytest.mark.ci
@pytest.mark.unit_test
def test_map_v6_particle_model_old_type_with_particle_group():
    """Old unit type with particle_type_000 subgroup should use v6 logic."""
    pt_group = _make_h5_group({
        'HAS_FILM_DIFFUSION': True,
        'HAS_PORE_DIFFUSION': True,
        'HAS_SURFACE_DIFFUSION': True,
    })
    group = _make_h5_group({'NPARTYPE': 1}, subgroups={'particle_type_000': pt_group})
    assert map_unit_to_particle_model('GENERAL_RATE_MODEL', group) == "1D (radial coordinate)"


# %% map_unit_type_to_column_model – v6 types

@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("unit_type, expected", [
    ('COLUMN_MODEL_1D', "1D (axial coordinate)"),
    ('COLUMN_MODEL_2D', "2D (axial and radial coordinate)"),
])
def test_map_unit_type_to_column_model_v6(unit_type, expected):
    """V6 unit types should map to correct column resolution."""
    assert map_unit_type_to_column_model(unit_type) == expected


# %% extract_config_data_from_unit – v6 paths

@pytest.mark.ci
@pytest.mark.unit_test
def test_extract_v6_config_grm_with_surface_diff_and_core():
    """V6 GRM-like config with surface diffusion and core radius > 0."""
    ads_group = _make_h5_group({'IS_KINETIC': True})
    pt_group = _make_h5_group(
        {
            'HAS_FILM_DIFFUSION': True,
            'HAS_PORE_DIFFUSION': True,
            'HAS_SURFACE_DIFFUSION': True,
            'ADSORPTION_MODEL': 'LINEAR',
            'NBOUND': np.array([1]),
            'PAR_CORERADIUS': 0.001,
        },
        subgroups={'adsorption': ads_group},
    )
    group = _make_h5_group(
        {
            'NPARTYPE': 1,
            'COL_POROSITY': 0.37,
            'COL_DISPERSION': 5.75e-08,
        },
        subgroups={'particle_type_000': pt_group},
    )

    config = extract_config_data_from_unit('COLUMN_MODEL_1D', group)

    assert config['PSD'] == "Yes"
    assert 'add_particles' not in config
    assert config['particle_resolution'] == "1D (radial coordinate)"
    assert config['particle_nonlimiting_filmDiff'] == "No"
    assert config['has_binding'] == "Yes"
    assert config['particle_has_surfDiff'] == "Yes"
    assert config['particle_has_core'] == "Yes"
    assert config['advanced_mode'] == "On"


@pytest.mark.ci
@pytest.mark.unit_test
def test_extract_v6_config_ptd_detection():
    """V6 config where binding model is an array with different entries should set PTD=Yes."""
    ads_group = _make_h5_group({'IS_KINETIC': True})
    pt_group = _make_h5_group(
        {
            'HAS_FILM_DIFFUSION': True,
            'HAS_PORE_DIFFUSION': False,
            'HAS_SURFACE_DIFFUSION': False,
            'ADSORPTION_MODEL': np.array([b'LINEAR', b'LANGMUIR']),
            'NBOUND': np.array([1]),
        },
        subgroups={'adsorption': ads_group},
    )
    group = _make_h5_group(
        {
            'NPARTYPE': 2,
            'COL_POROSITY': 0.37,
            'COL_DISPERSION': 5.75e-08,
        },
        subgroups={'particle_type_000': pt_group},
    )

    config = extract_config_data_from_unit('COLUMN_MODEL_1D', group)

    assert config['PSD'] == "Particle size distribution"
    assert config['advanced_mode'] == "On"


@pytest.mark.ci
@pytest.mark.unit_test
def test_map_unit_type_to_column_model_3d():
    """A 3D unit type should map to 3D resolution."""
    assert map_unit_type_to_column_model('GENERAL_RATE_MODEL_3D') == "3D (axial, radial and angular coordinate)"


@pytest.mark.ci
@pytest.mark.unit_test
def test_extract_config_0d_with_flow_filter():
    """CSTR with FLOWRATE_FILTER > 0 should set has_filter to Yes."""
    group = _make_h5_group({
        'UNIT_TYPE': b'CSTR',
        'NPARTYPE': 0,
        'TOTAL_POROSITY': 1.0,
        'FLOWRATE_FILTER': 1.0,
    })
    config = extract_config_data_from_unit('CSTR', group)
    assert config['has_filter'] == "Yes"


@pytest.mark.ci
@pytest.mark.unit_test
def test_extract_config_nbound_greater_than_one():
    """NBOUND > 1 should set has_mult_bnd_states to Yes."""
    ads_group = _make_h5_group({'IS_KINETIC': True})
    pt_group = _make_h5_group(
        {
            'HAS_FILM_DIFFUSION': True,
            'HAS_PORE_DIFFUSION': True,
            'HAS_SURFACE_DIFFUSION': False,
            'ADSORPTION_MODEL': b'LINEAR',
            'NBOUND': 2,
            'PAR_CORERADIUS': 0.001,
        },
        subgroups={'adsorption': ads_group},
    )
    group = _make_h5_group(
        {
            'NPARTYPE': 1,
            'COL_POROSITY': 0.37,
            'COL_DISPERSION': 5.75e-08,
        },
        subgroups={'particle_type_000': pt_group},
    )
    config = extract_config_data_from_unit('COLUMN_MODEL_1D', group)
    assert config['has_mult_bnd_states'] == "Yes"


@pytest.mark.ci
@pytest.mark.unit_test
def test_map_unit_to_particle_model_column_model_2d_no_particle():
    """COLUMN_MODEL_2D without particles (NPARTYPE=0) should return None."""
    group = _make_h5_group({'NPARTYPE': 0})
    assert map_unit_to_particle_model('COLUMN_MODEL_2D', group) is None


@pytest.mark.ci
@pytest.mark.unit_test
def test_extract_v6_config_psd_shared_surfDiff():
    """V6 PSD config should use shared surfDiff key (not per-type) since PSD uses shared config."""
    ads_group = _make_h5_group({'IS_KINETIC': True})
    pt0_group = _make_h5_group(
        {
            'HAS_FILM_DIFFUSION': True,
            'HAS_PORE_DIFFUSION': True,
            'HAS_SURFACE_DIFFUSION': True,
            'ADSORPTION_MODEL': b'LINEAR',
            'NBOUND': np.array([1]),
        },
        subgroups={'adsorption': ads_group},
    )
    group = _make_h5_group(
        {
            'NPARTYPE': 2,
            'COL_POROSITY': 0.37,
            'COL_DISPERSION': 5.75e-08,
        },
        subgroups={'particle_type_000': pt0_group},
    )

    config = extract_config_data_from_unit('COLUMN_MODEL_1D', group)

    assert config['PSD'] == "Particle size distribution"
    assert config['particle_has_surfDiff'] == "Yes"
    assert 'parType_1_has_surfDiff' not in config  # shared config, not per-type


@pytest.mark.ci
@pytest.mark.unit_test
def test_extract_v6_config_single_particle_no_binding():
    """V6 single particle without binding should set has_binding to No."""
    pt_group = _make_h5_group(
        {
            'HAS_FILM_DIFFUSION': True,
            'HAS_PORE_DIFFUSION': True,
            'HAS_SURFACE_DIFFUSION': False,
            'ADSORPTION_MODEL': b'NONE',
        },
        subgroups={},
    )
    group = _make_h5_group(
        {
            'NPARTYPE': 1,
            'COL_POROSITY': 0.37,
            'COL_DISPERSION': 5.75e-08,
        },
        subgroups={'particle_type_000': pt_group},
    )

    config = extract_config_data_from_unit('COLUMN_MODEL_1D', group)

    assert config['add_particles'] == "Yes"
    assert config['has_binding'] == "No"
    assert 'particle_has_surfDiff' not in config


# %% extract_config_data_from_unit – column geometry

@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("unit_type, expected_geometry", [
    ('GENERAL_RATE_MODEL',            "Axial flow cylinder"),
    ('COLUMN_MODEL_1D',               "Axial flow cylinder"),
    ('CSTR',                          "Axial flow cylinder"),
    ('RADIAL_GENERAL_RATE_MODEL',     "Radial flow cylinder"),
    ('RADIAL_COLUMN_MODEL_1D',        "Radial flow cylinder"),
    ('FRUSTUM_GENERAL_RATE_MODEL',    "Frustum"),
    ('FRUSTUM_COLUMN_MODEL_1D',       "Frustum"),
])
def test_extract_config_sets_column_type(unit_type, expected_geometry):
    """extract_config_data_from_unit should set column_type from the CADET unit type."""
    group = _make_h5_group({
        'COL_POROSITY': 1.0,
        'TOTAL_POROSITY': 1.0,
    })
    config = extract_config_data_from_unit(unit_type, group)
    assert config['column_type'] == expected_geometry


# %% Binding model extraction from h5

@pytest.mark.ci
@pytest.mark.unit_test
@pytest.mark.parametrize("cadet_name, expected", [
    ('LINEAR', 'Linear'),
    ('MULTI_COMPONENT_LANGMUIR', 'Langmuir'),
    ('STERIC_MASS_ACTION', 'SMA'),
])
def test_extract_binding_model_v6_known_models(cadet_name, expected):
    """Known CADET binding model names should map to the generator's model name."""
    ads_group = _make_h5_group({'IS_KINETIC': True})
    pt_group = _make_h5_group(
        {
            'HAS_FILM_DIFFUSION': True,
            'HAS_PORE_DIFFUSION': True,
            'HAS_SURFACE_DIFFUSION': False,
            'ADSORPTION_MODEL': cadet_name.encode(),
            'NBOUND': np.array([1]),
        },
        subgroups={'adsorption': ads_group},
    )
    group = _make_h5_group(
        {
            'NPARTYPE': 1,
            'COL_POROSITY': 0.37,
            'COL_DISPERSION': 5.75e-08,
        },
        subgroups={'particle_type_000': pt_group},
    )

    config = extract_config_data_from_unit('COLUMN_MODEL_1D', group)

    assert config['binding_model'] == expected


@pytest.mark.ci
@pytest.mark.unit_test
def test_extract_binding_model_v6_unknown_falls_back_to_arbitrary():
    """An unknown CADET binding model should fall back to 'Arbitrary'."""
    ads_group = _make_h5_group({'IS_KINETIC': False})
    pt_group = _make_h5_group(
        {
            'HAS_FILM_DIFFUSION': True,
            'HAS_PORE_DIFFUSION': False,
            'HAS_SURFACE_DIFFUSION': False,
            'ADSORPTION_MODEL': b'MOBILE_PHASE_MODULATOR',
            'NBOUND': np.array([1]),
        },
        subgroups={'adsorption': ads_group},
    )
    group = _make_h5_group(
        {
            'NPARTYPE': 1,
            'COL_POROSITY': 0.37,
            'COL_DISPERSION': 5.75e-08,
        },
        subgroups={'particle_type_000': pt_group},
    )

    config = extract_config_data_from_unit('COLUMN_MODEL_1D', group)

    assert config['binding_model'] == "Arbitrary"


@pytest.mark.ci
@pytest.mark.unit_test
def test_extract_binding_model_v5_linear():
    """V5 GRM with LINEAR binding should extract binding_model='Linear'."""
    ads_group = _make_h5_group({'IS_KINETIC': True})
    group = _make_h5_group(
        {
            'COL_POROSITY': 0.37,
            'COL_DISPERSION': 5.75e-08,
            'ADSORPTION_MODEL': b'LINEAR',
            'NBOUND': np.array([1]),
            'PAR_SURFDIFFUSION': 1.0e-10,
        },
        subgroups={'adsorption': ads_group},
    )

    config = extract_config_data_from_unit('GENERAL_RATE_MODEL', group)

    assert config['binding_model'] == "Linear"


@pytest.mark.ci
@pytest.mark.unit_test
def test_extract_binding_model_not_set_when_no_binding():
    """When ADSORPTION_MODEL is NONE, binding_model should not be in config."""
    pt_group = _make_h5_group(
        {
            'HAS_FILM_DIFFUSION': True,
            'HAS_PORE_DIFFUSION': True,
            'HAS_SURFACE_DIFFUSION': False,
            'ADSORPTION_MODEL': b'NONE',
        },
        subgroups={},
    )
    group = _make_h5_group(
        {
            'NPARTYPE': 1,
            'COL_POROSITY': 0.37,
            'COL_DISPERSION': 5.75e-08,
        },
        subgroups={'particle_type_000': pt_group},
    )

    config = extract_config_data_from_unit('COLUMN_MODEL_1D', group)

    assert 'binding_model' not in config


@pytest.mark.ci
@pytest.mark.unit_test
def test_extract_binding_model_v6_multi_partype_uses_first():
    """V6 multi-particle config should use first particle type's binding model as shared config."""
    ads_group_0 = _make_h5_group({'IS_KINETIC': True})
    pt0_group = _make_h5_group(
        {
            'HAS_FILM_DIFFUSION': True,
            'HAS_PORE_DIFFUSION': False,
            'HAS_SURFACE_DIFFUSION': False,
            'ADSORPTION_MODEL': b'STERIC_MASS_ACTION',
            'NBOUND': np.array([1]),
        },
        subgroups={'adsorption': ads_group_0},
    )
    pt1_group = _make_h5_group(
        {
            'HAS_FILM_DIFFUSION': True,
            'HAS_PORE_DIFFUSION': False,
            'HAS_SURFACE_DIFFUSION': False,
            'ADSORPTION_MODEL': b'LINEAR',
            'NBOUND': np.array([1]),
        },
        subgroups={'adsorption': ads_group_0},
    )
    group = _make_h5_group(
        {
            'NPARTYPE': 2,
            'COL_POROSITY': 0.37,
            'COL_DISPERSION': 5.75e-08,
        },
        subgroups={
            'particle_type_000': pt0_group,
            'particle_type_001': pt1_group,
        },
    )

    config = extract_config_data_from_unit('COLUMN_MODEL_1D', group)

    # Shared binding_model comes from first particle type
    assert config['binding_model'] == "SMA"
    # Per-particle-type keys are not extracted (only used in dev mode)
    assert 'parType_1_binding_model' not in config
    assert 'parType_2_binding_model' not in config


# %% Crystallization detection and extraction

@pytest.mark.ci
@pytest.mark.unit_test
def test_is_crystallization_unit_true():
    group = _make_h5_group({'REACTION_MODEL': b'CRYSTALLIZATION'})
    assert _is_crystallization_unit(group) is True


@pytest.mark.ci
@pytest.mark.unit_test
def test_is_crystallization_unit_false():
    group = _make_h5_group({'REACTION_MODEL': b'NONE'})
    assert _is_crystallization_unit(group) is False


@pytest.mark.ci
@pytest.mark.unit_test
def test_is_crystallization_unit_missing():
    group = _make_h5_group({})
    assert _is_crystallization_unit(group) is False


@pytest.mark.ci
@pytest.mark.unit_test
def test_get_crystallization_reaction_group_reaction_bulk():
    reaction_bulk = _make_h5_group({'CRY_MODE': 1})
    group = _make_h5_group({}, subgroups={'reaction_bulk': reaction_bulk})
    assert _get_crystallization_reaction_group(group) is reaction_bulk


@pytest.mark.ci
@pytest.mark.unit_test
def test_get_crystallization_reaction_group_reaction():
    reaction = _make_h5_group({'CRY_MODE': 3})
    group = _make_h5_group({}, subgroups={'reaction': reaction})
    assert _get_crystallization_reaction_group(group) is reaction


@pytest.mark.ci
@pytest.mark.unit_test
def test_get_crystallization_reaction_group_none():
    group = _make_h5_group({})
    assert _get_crystallization_reaction_group(group) is None


@pytest.mark.ci
@pytest.mark.unit_test
def test_extract_crystallization_config_no_reaction_group():
    group = _make_h5_group({'REACTION_MODEL': b'CRYSTALLIZATION'})
    config = extract_crystallization_config('CSTR', group)
    assert config['model_type'] == 'Crystallization'
    assert config['cry_column_type'] == 'CSTR'
    assert config['cry_has_primary_formation'] == 'Yes'
    assert config['cry_has_aggregation'] == 'No'
    assert config['cry_has_fragmentation'] == 'No'


@pytest.mark.ci
@pytest.mark.unit_test
def test_extract_crystallization_config_cstr_primary():
    reaction = _make_h5_group({
        'CRY_GROWTH_RATE_CONSTANT': 5e-6,
        'CRY_GROWTH_DISPERSION_RATE': 2e-14,
        'CRY_PRIMARY_NUCLEATION_RATE': 1e6,
        'CRY_SECONDARY_NUCLEATION_RATE': 0.0,
        'CRY_P': 0.0,
    })
    group = _make_h5_group(
        {'REACTION_MODEL': b'CRYSTALLIZATION'},
        subgroups={'reaction_bulk': reaction},
    )
    config = extract_crystallization_config('CSTR', group)
    assert config['cry_column_type'] == 'CSTR'
    assert config['cry_has_primary_formation'] == 'Yes'
    assert config['cry_has_growth_dispersion'] == 'Yes'
    assert config['cry_has_secondary_nucleation'] == 'No'
    assert config['cry_size_dependent_growth'] == 'No'
    assert config['cry_has_aggregation'] == 'No'
    assert config['cry_has_fragmentation'] == 'No'


@pytest.mark.ci
@pytest.mark.unit_test
def test_extract_crystallization_config_dpfr_aggregation():
    reaction = _make_h5_group({
        'CRY_MODE': 3,
        'CRY_GROWTH_RATE_CONSTANT': 5e-6,
        'CRY_GROWTH_DISPERSION_RATE': 2.5e-15,
        'CRY_PRIMARY_NUCLEATION_RATE': 5.0,
        'CRY_SECONDARY_NUCLEATION_RATE': 4e8,
        'CRY_P': 0.0,
        'CRY_AGGREGATION_INDEX': 0,
    })
    group = _make_h5_group(
        {'REACTION_MODEL': b'CRYSTALLIZATION', 'COL_DISPERSION': 4.2e-5},
        subgroups={'reaction': reaction},
    )
    config = extract_crystallization_config('LUMPED_RATE_MODEL_WITHOUT_PORES', group)
    assert config['cry_column_type'] == 'DPFR'
    assert config['cry_has_primary_formation'] == 'Yes'
    assert config['cry_has_aggregation'] == 'Yes'
    assert config['cry_aggregation_kernel'] == 'Constant'
    assert config['cry_has_fragmentation'] == 'No'
    assert config['cry_has_axial_dispersion'] == 'Yes'


@pytest.mark.ci
@pytest.mark.unit_test
def test_extract_crystallization_config_agg_without_index():
    reaction = _make_h5_group({'CRY_MODE': 3})
    group = _make_h5_group(
        {'REACTION_MODEL': b'CRYSTALLIZATION'},
        subgroups={'reaction_bulk': reaction},
    )
    config = extract_crystallization_config('CSTR', group)
    assert config['cry_has_aggregation'] == 'Yes'
    assert config['cry_aggregation_kernel'] == 'Constant'


@pytest.mark.ci
@pytest.mark.unit_test
def test_extract_crystallization_config_dpfr_no_axial_dispersion():
    reaction = _make_h5_group({})
    group = _make_h5_group(
        {'REACTION_MODEL': b'CRYSTALLIZATION', 'COL_DISPERSION': 0.0},
        subgroups={'reaction': reaction},
    )
    config = extract_crystallization_config('LUMPED_RATE_MODEL_WITHOUT_PORES', group)
    assert config['cry_has_axial_dispersion'] == 'No'


# %% Reaction config extraction tests

@pytest.mark.ci
@pytest.mark.unit_test
def test_extract_reaction_config_from_bulk_group():
    """Test _extract_reaction_config reading from reaction_bulk subgroup."""
    reaction_bulk = _make_h5_group({'REACTION_MODEL': b'MASS_ACTION_LAW'})
    unit_group = _make_h5_group({}, subgroups={'reaction_bulk': reaction_bulk})
    config = {}
    _extract_reaction_config(config, unit_group)
    assert config['reaction_model'] == 'Mass Action Law'
    assert config['has_reaction_bulk'] == "Yes"


@pytest.mark.ci
@pytest.mark.unit_test
def test_extract_reaction_config_from_unit_level():
    """Test _extract_reaction_config reading REACTION_MODEL at unit level."""
    unit_group = _make_h5_group({'REACTION_MODEL': b'MICHAELIS_MENTEN'})
    config = {}
    _extract_reaction_config(config, unit_group)
    assert config['reaction_model'] == 'Michaelis Menten'
    assert config['has_reaction_bulk'] == "Yes"


@pytest.mark.ci
@pytest.mark.unit_test
def test_extract_reaction_config_unknown_model():
    """Test _extract_reaction_config with unknown model maps to Arbitrary."""
    unit_group = _make_h5_group({'REACTION_MODEL': b'SOME_UNKNOWN_MODEL'})
    config = {}
    _extract_reaction_config(config, unit_group)
    assert config['reaction_model'] == 'Arbitrary'
    assert config['has_reaction_bulk'] == "Yes"


@pytest.mark.ci
@pytest.mark.unit_test
def test_extract_reaction_config_none_model():
    """Test _extract_reaction_config with NONE does not set reaction config."""
    unit_group = _make_h5_group({'REACTION_MODEL': b'NONE'})
    config = {}
    _extract_reaction_config(config, unit_group)
    assert 'reaction_model' not in config
    assert 'has_reaction_bulk' not in config
