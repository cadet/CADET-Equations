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
    map_unit_type_to_column_model,
    map_unit_to_particle_model,
    CADET_column_unit_types,
)


# %% Helpers

def _make_h5_group(mapping: dict):
    """Create a mock HDF5 group that mimics h5py dataset access.

    Keys present in *mapping* return mock datasets whose [()] call
    yields the corresponding value.  Missing keys return None.
    """
    group = MagicMock()

    def get_side_effect(k):
        if k in mapping:
            ds = MagicMock()
            val = mapping[k]
            type(ds).__getitem__ = lambda self, idx: val
            return ds
        return None

    group.get = get_side_effect
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
    ('GENERAL_RATE_MODEL',                 "1D (axial coordinate)"),
    ('GENERAL_RATE_MODEL_DG',              "1D (axial coordinate)"),
    ('LUMPED_RATE_MODEL_WITHOUT_PORES',    "1D (axial coordinate)"),
    ('LUMPED_RATE_MODEL_WITH_PORES',       "1D (axial coordinate)"),
    ('LUMPED_RATE_MODEL_WITHOUT_PORES_DG', "1D (axial coordinate)"),
    ('LUMPED_RATE_MODEL_WITH_PORES_DG',    "1D (axial coordinate)"),
    ('GENERAL_RATE_MODEL_2D',              "2D (axial and radial coordinate)"),
    ('CSTR',                               "0D (Homogeneous Tank)"),
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
    """The constant should contain all 8 supported CADET column unit types."""
    expected = [
        'GENERAL_RATE_MODEL', 'LUMPED_RATE_MODEL_WITHOUT_PORES',
        'LUMPED_RATE_MODEL_WITH_PORES', 'GENERAL_RATE_MODEL_DG',
        'LUMPED_RATE_MODEL_WITHOUT_PORES_DG', 'LUMPED_RATE_MODEL_WITH_PORES_DG',
        'GENERAL_RATE_MODEL_2D', 'CSTR'
    ]
    assert len(CADET_column_unit_types) == 8
    for t in expected:
        assert t in CADET_column_unit_types
