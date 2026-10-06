"""
Helpers to extract a generator configuration from CADET HDF5 files.

This module provides minimal mapping functions used when a user
uploads a CADET HDF5 file to pre-populate the generator UI.
"""

import re

import h5py
import streamlit as st

CADET_binding_model_map = {
    "LINEAR": "Linear",
    "MULTI_COMPONENT_LANGMUIR": "Langmuir",
    "STERIC_MASS_ACTION": "SMA",
}

CADET_reaction_model_map = {
    "MASS_ACTION_LAW": "Mass Action Law",
    "MASS_ACTION_LAW_CROSS_PHASE": "Mass Action Law",
    "MICHAELIS_MENTEN": "Michaelis Menten",
}

CADET_particle_geometry_map = {
    "SPHERE": "Sphere",
    "CYLINDER": "Cylinder",
    "SLAB": "Slab",
}

# Binding widgets are not prefixed with particle_ when a single configuration is shared.
_BINDING_KEYS = {"binding_model", "req_binding", "has_mult_bnd_states"}

# Marks a configuration that cannot be expressed outside developer mode.
_DEV_MODE_REQUIRED = "_dev_mode_required"

# Session state key of the particle type count, which developer mode uses in place of
# the "Add particles" selectbox.
_N_PAR_TYPE_KEY = "N^\\mathrm{p}"


# Column geometry is a dedicated field since CADET-Core v6; before that it was encoded
# in the unit type name. Keys are the values of GEOMETRY in /input/model/unit_XXX.
CADET_geometry_map = {
    "AXIAL_FLOW_CYLINDER": "Axial flow cylinder",
    "RADIAL_FLOW_CYLINDER_SHELL": "Radial flow cylinder",
    "AXIAL_FLOW_FRUSTUM": "Frustum",
}

# Geometries CADET-Core supports that have no counterpart in CADET-Equations.
CADET_unsupported_geometries = {"SMOOTHLY_VARYING"}

CADET_column_unit_types = [
    "CSTR",
    "COLUMN_MODEL_1D",
    "COLUMN_MODEL_2D",
    "GENERAL_RATE_MODEL",
    "GRM",
    "LUMPED_RATE_MODEL_WITHOUT_PORES",
    "LRM",
    "LUMPED_RATE_MODEL_WITH_PORES",
    "LRMP",
    "DISPERSIVE_PLUG_FLOW_REACTOR",
    "DPFR",
    "GENERAL_RATE_MODEL_2D",
    # unit operation identifiers, which encode discretization and (legacy) geometry
    "AXIAL_COLUMN_MODEL_1D_COLLOCATION_DG",
    "AXIAL_COLUMN_MODEL_1D_FV",
    "RADIAL_COLUMN_MODEL_1D_FV",
    "FRUSTUM_COLUMN_MODEL_1D_FV",
    "VARIABLE_CROSS_SECTION_COLUMN_MODEL_1D_DG",
    "GENERAL_RATE_MODEL_FV",
    "RADIAL_GENERAL_RATE_MODEL_FV",
    "FRUSTUM_GENERAL_RATE_MODEL_FV",
    "LUMPED_RATE_MODEL_WITH_PORES_FV",
    "RADIAL_LUMPED_RATE_MODEL_WITH_PORES_FV",
    "FRUSTUM_LUMPED_RATE_MODEL_WITH_PORES_FV",
    "LUMPED_RATE_MODEL_WITHOUT_PORES_FV",
    "RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES_FV",
    "FRUSTUM_LUMPED_RATE_MODEL_WITHOUT_PORES_FV",
    "LUMPED_RATE_MODEL_WITHOUT_PORES_COLLOCATIONDG",
    # old interface
    "FRUSTUM_COLUMN_MODEL_1D",
    "FRUSTUM_GENERAL_RATE_MODEL",
    "FRUSTUM_LUMPED_RATE_MODEL_WITHOUT_PORES",
    "FRUSTUM_LUMPED_RATE_MODEL_WITH_PORES",
    "RADIAL_COLUMN_MODEL_1D",
    "RADIAL_GENERAL_RATE_MODEL",
    "RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES",
    "RADIAL_LUMPED_RATE_MODEL_WITH_PORES",
    "GENERAL_RATE_MODEL_DG",
    "LUMPED_RATE_MODEL_WITHOUT_PORES_DG",
    "LUMPED_RATE_MODEL_WITH_PORES_DG",
]


def is_v6_interface(unit_type, h5_unit_group):
    """Detect a v6-style CADET HDF5 unit group.

    Returns True when the file structure matches the newer v6 layout.
    """
    if re.search("COLUMN_MODEL", unit_type):
        return True
    nParType = get_h5_value(h5_unit_group, "NPARTYPE")
    return bool(nParType is not None and nParType >= 1 and "particle_type_000" in h5_unit_group)


def get_h5_value(unit_group, key: str, firstEntryIfList=True):
    """Read a value from an HDF5 group and return a native Python type.

    The function decodes bytes and returns the first element of array-like
    values when appropriate.
    """

    value = unit_group.get(key)

    if value is None:
        return None

    value = value[()]

    if isinstance(value, bytes):
        value = value.decode("utf-8")

    # If the value is a NumPy array or similar, get the first element
    if hasattr(value, "__len__") and not isinstance(value, (str, bytes)):
        if len(value) > 0 and firstEntryIfList:
            value = value[0]

    return value


def map_unit_type_to_column_geometry(cadet_unit_type, h5_unit_group=None):
    """Return the generator's column geometry label for a CADET unit.

    Since CADET-Core v6 the geometry is a dedicated GEOMETRY field and the unit type
    no longer carries it, so GEOMETRY takes precedence over the unit type name.
    """

    if re.search("CSTR", cadet_unit_type):
        return "Mixed tank"

    geometry = get_h5_value(h5_unit_group, "GEOMETRY") if h5_unit_group is not None else None

    if geometry is not None:
        if geometry in CADET_geometry_map:
            return CADET_geometry_map[geometry]
        if geometry in CADET_unsupported_geometries:
            st.sidebar.warning(
                f"Column geometry {geometry} is not available in CADET-Equations, defaulting to an axial flow cylinder."
            )
            return "Axial flow cylinder"
        st.sidebar.warning(f"Unknown column geometry {geometry}, defaulting to an axial flow cylinder.")
        return "Axial flow cylinder"

    # Pre-v6 files encode the geometry in the unit type name
    if re.search("RADIAL", cadet_unit_type):
        return "Radial flow cylinder"
    elif re.search("FRUSTUM", cadet_unit_type):
        return "Frustum"
    else:
        return "Axial flow cylinder"


def map_unit_type_to_column_model(cadet_unit_type, h5_unit_group=None):
    """Map a CADET unit to the generator's column resolution label."""

    if re.search("3D", cadet_unit_type):
        return "3D (axial, radial and angular coordinate)"
    elif re.search("2D", cadet_unit_type):
        return "2D (axial and radial coordinate)"
    elif re.search("CSTR", cadet_unit_type):
        return "0D (Homogeneous Tank)"
    elif cadet_unit_type not in CADET_column_unit_types:
        raise ValueError(f"Invalid unit type: {cadet_unit_type}. Must be one of {CADET_column_unit_types}.")

    # A radial flow column resolves the radial rather than the axial coordinate
    if map_unit_type_to_column_geometry(cadet_unit_type, h5_unit_group) == "Radial flow cylinder":
        return "1D (radial coordinate)"

    return "1D (axial coordinate)"


def map_unit_to_particle_model(cadet_unit_type, h5_unit_group):
    """Return a textual particle model mapping for a CADET unit type/group."""

    if is_v6_interface(cadet_unit_type, h5_unit_group):
        return _map_v6_particle_model(h5_unit_group)

    if re.search("WITHOUT_PORES", cadet_unit_type) or re.search("CSTR", cadet_unit_type):
        if get_h5_value(h5_unit_group, "TOTAL_POROSITY") == 1.0:
            return None
        if get_h5_value(h5_unit_group, "CONST_SOLID_VOLUME") == 0.0:
            return None

    elif get_h5_value(h5_unit_group, "COL_POROSITY") == 1.0:
        return None

    if re.search("GENERAL_RATE", cadet_unit_type):
        return "1D (radial coordinate)"
    elif re.search("LUMPED_RATE", cadet_unit_type) or re.search("CSTR", cadet_unit_type):
        return "0D (homogeneous)"
    elif cadet_unit_type in CADET_column_unit_types:
        return None
    else:
        raise ValueError(f"Invalid unit type: {cadet_unit_type}. Must be one of {CADET_column_unit_types}.")


def _particle_resolution(pt_group):
    """Return the generator's particle resolution for one particle type group.

    Mirrors the particle transport type CADET-Core derives from the HAS_* flags: only a
    general rate particle resolves the particle radius, the homogeneous and equilibrium
    particles do not.
    """

    has_pore_diff = get_h5_value(pt_group, "HAS_PORE_DIFFUSION")
    has_surf_diff = get_h5_value(pt_group, "HAS_SURFACE_DIFFUSION")

    if has_pore_diff or has_surf_diff:
        return "1D (radial coordinate)"

    return "0D (homogeneous)"


def _map_v6_particle_model(h5_unit_group):
    """Determine particle model for v6 interface using HAS_* flags in particle_type_000."""

    nParType = get_h5_value(h5_unit_group, "NPARTYPE")
    if nParType is None or nParType < 1:
        return None

    pt_group = h5_unit_group.get("particle_type_000")
    if pt_group is None:
        return None

    return _particle_resolution(pt_group)


def extract_config_data_from_unit(unit_type, h5_unit_group):

    config = {}

    config["advanced_mode"] = "Off"
    config["var_format"] = "CADET"
    config["show_eq_description"] = True
    config["model_assumptions"] = True

    config["column_type"] = map_unit_type_to_column_geometry(unit_type, h5_unit_group)
    config["column_resolution"] = map_unit_type_to_column_model(unit_type, h5_unit_group)

    if re.search("0D", config["column_resolution"]):
        flow_filter = get_h5_value(h5_unit_group, "FLOWRATE_FILTER")
        config["has_filter"] = "No"
        if flow_filter is not None:
            config["has_filter"] = "Yes" if flow_filter > 0.0 else "No"

    if re.search("2D", config["column_resolution"]):
        # The 2D models renamed COL_DISPERSION to COL_DISPERSION_AXIAL
        Dax = get_h5_value(h5_unit_group, "COL_DISPERSION_AXIAL")
        if Dax is None:
            Dax = get_h5_value(h5_unit_group, "COL_DISPERSION")
        if Dax is not None:
            config["has_axial_dispersion"] = "No" if Dax < 1e-20 else "Yes"

        Drad = get_h5_value(h5_unit_group, "COL_DISPERSION_RADIAL")
        if Drad is not None:
            config["has_radial_dispersion"] = "No" if Drad < 1e-20 else "Yes"

    elif re.search("1D", config["column_resolution"]):
        Dax = get_h5_value(h5_unit_group, "COL_DISPERSION")
        if Dax is not None:
            config["has_axial_dispersion"] = "No" if Dax < 1e-20 else "Yes"

    par_model = map_unit_to_particle_model(unit_type, h5_unit_group)

    if par_model is not None:
        config["add_particles"] = "Yes"

        config["particle_resolution"] = par_model

        if is_v6_interface(unit_type, h5_unit_group):
            _extract_v6_particle_config(config, h5_unit_group, par_model)
        else:
            _extract_v5_particle_config(config, unit_type, h5_unit_group, par_model)

    else:
        config["add_particles"] = "No"

    dev_mode_required = config.pop(_DEV_MODE_REQUIRED, False)

    if dev_mode_required:
        # Developer mode implies the advanced options and replaces both the "Add particles"
        # and the particle size distribution selectbox with a particle type count
        config["advanced_mode"] = "On"
        if _N_PAR_TYPE_KEY not in config:
            n_par_type = get_h5_value(h5_unit_group, "NPARTYPE")
            if n_par_type is None:
                n_par_type = 1 if config.get("add_particles") == "Yes" else 0
            config[_N_PAR_TYPE_KEY] = int(n_par_type)
        config.pop("add_particles", None)
        config.pop("PSD", None)

    elif config["advanced_mode"] == "On":
        # Advanced mode uses single PSD selectbox with 3 options
        add_par = config.pop("add_particles", "No")
        if add_par == "Yes":
            if config.get("PSD") == "Yes":
                config["PSD"] = "Particle size distribution"
            else:
                config["PSD"] = "Yes"
        else:
            config["PSD"] = "No"
    else:
        config.pop("dev_mode", None)
        config.pop("PSD", None)
        config.pop("particle_has_core", None)
        config.pop("has_radial_dispersion", None)
        config.pop("has_mult_bnd_states", None)
        config.pop("particle_geometry", None)

    return config


def _extract_v5_particle_config(config, unit_type, h5_unit_group, par_model):
    """Extract particle configuration from v5 interface (particle info at unit level)."""

    nonlimiting = bool(re.search("WITHOUT_PORES", unit_type))

    nParType = get_h5_value(h5_unit_group, "NPARTYPE")
    nParType = 1 if nParType is None else nParType

    if nParType > 1:
        config["advanced_mode"] = "On"
        config["PSD"] = "Yes"

    config["particle_nonlimiting_filmDiff"] = "Yes" if nonlimiting else "No"

    binding_model = get_h5_value(h5_unit_group, "ADSORPTION_MODEL", firstEntryIfList=False)

    config["has_binding"] = "No"

    if binding_model is not None:
        if not isinstance(binding_model, str):
            if len(binding_model) > 1:
                binding_model = binding_model[0]

        if binding_model != "NONE":
            config["has_binding"] = "Yes"

            config["binding_model"] = CADET_binding_model_map.get(binding_model, "Arbitrary")
            if binding_model not in CADET_binding_model_map:
                st.sidebar.warning(
                    f"Binding model {binding_model} not implemented in CADET-Equations, default to arbitrary binding"
                )

            if nParType > 1:
                ads_group = h5_unit_group["adsorption_000"]
            else:
                ads_group = h5_unit_group["adsorption"]

            config["req_binding"] = "Kinetic" if get_h5_value(ads_group, "IS_KINETIC") else "Rapid-equilibrium"
            if config["binding_model"] == "Arbitrary":
                config["has_mult_bnd_states"] = "No"
                if get_h5_value(ads_group, "NBOUND") is not None:
                    config["has_mult_bnd_states"] = "Yes" if get_h5_value(ads_group, "NBOUND") > 1 else "No"

            if par_model == "1D (radial coordinate)":
                surfDiff = get_h5_value(h5_unit_group, "PAR_SURFDIFFUSION")
                config["particle_has_surfDiff"] = "Yes" if surfDiff is not None and surfDiff > 0.0 else "No"

    if par_model == "1D (radial coordinate)":
        _extract_particle_core_config(config, h5_unit_group)

    _extract_reaction_config(config, h5_unit_group)


def particle_type_groups(h5_unit_group, n_par_type):
    """Return the particle_type_XXX groups of a unit, in index order.

    Stops at the first missing group, so a file that announces more types than it
    stores yields only the ones actually present.
    """

    groups = []
    for j in range(n_par_type):
        group = h5_unit_group.get(f"particle_type_{j:03d}")
        if group is None:
            break
        groups.append(group)
    return groups


def _extract_v6_particle_config(config, h5_unit_group, par_model):
    """Extract particle configuration from v6 interface (particle info in particle_type_xxx subgroups).

    Particle types that share all settings CADET-Equations models are a particle size
    distribution and keep the shared widgets. Types that genuinely differ are written to
    their own parType_X_ keys, which the generator only offers in developer mode.
    """

    nParType = get_h5_value(h5_unit_group, "NPARTYPE")
    nParType = 1 if nParType is None else int(nParType)

    pt_groups = particle_type_groups(h5_unit_group, nParType)

    if not pt_groups:
        return

    if len(pt_groups) < nParType:
        st.sidebar.warning(
            f"The file announces {nParType} particle types but only stores {len(pt_groups)}; "
            "the stored configuration is applied to all of them."
        )

    type_configs = [_particle_type_config(pt_group, h5_unit_group) for pt_group in pt_groups]

    several_types = nParType > 1
    types_differ = any(type_config != type_configs[0] for type_config in type_configs[1:])

    # Flags describe the model as a whole rather than a single particle type
    for type_config in type_configs:
        for flag in ("dev_mode", _DEV_MODE_REQUIRED, "advanced_mode"):
            if flag in type_config:
                config[flag] = type_config.pop(flag)

    binding = [type_config.pop("has_binding", "No") for type_config in type_configs]
    config["has_binding"] = "Yes" if "Yes" in binding else "No"

    if types_differ:
        # distinct particle types can only be configured in developer mode
        config["dev_mode"] = True
        config[_DEV_MODE_REQUIRED] = True
        config[_N_PAR_TYPE_KEY] = len(type_configs)
        # the shared particle widgets are replaced by the per-type ones
        config.pop("particle_resolution", None)
        for j, type_config in enumerate(type_configs):
            for key, value in type_config.items():
                config[f"parType_{j + 1}_{key}"] = value
    else:
        if several_types:
            # types that differ only in size are a particle size distribution
            config["advanced_mode"] = "On"
            config["PSD"] = "Yes"
        for key, value in type_configs[0].items():
            config[key if key in _BINDING_KEYS else "particle_" + key] = value

    # Reactions are configured once; particle reactions are read from the first type
    _extract_reaction_config(config, h5_unit_group)
    _extract_particle_reaction_config(config, pt_groups[0])


def _particle_type_config(pt_group, h5_unit_group):
    """Return the generator settings of one particle type, keyed without any prefix."""

    config = {}
    resolution = _particle_resolution(pt_group)
    config["resolution"] = resolution

    has_film_diff = get_h5_value(pt_group, "HAS_FILM_DIFFUSION")
    config["nonlimiting_filmDiff"] = "No" if has_film_diff else "Yes"

    binding_model = get_h5_value(pt_group, "ADSORPTION_MODEL", firstEntryIfList=False)

    if binding_model is not None:
        if not isinstance(binding_model, str):
            if len(binding_model) > 1:
                binding_model = binding_model[0]

        if binding_model != "NONE":
            config["has_binding"] = "Yes"

            mapped_binding = CADET_binding_model_map.get(binding_model, "Arbitrary")
            config["binding_model"] = mapped_binding
            if binding_model not in CADET_binding_model_map:
                st.sidebar.warning(
                    f"Binding model {binding_model} not implemented in CADET-Equations, default to arbitrary binding"
                )

            ads_group = pt_group["adsorption"]
            config["req_binding"] = "Kinetic" if get_h5_value(ads_group, "IS_KINETIC") else "Rapid-equilibrium"
            if mapped_binding == "Arbitrary":
                nbound = get_h5_value(pt_group, "NBOUND")
                config["has_mult_bnd_states"] = "Yes" if nbound is not None and nbound > 1 else "No"

            if resolution == "1D (radial coordinate)":
                has_surf_diff = get_h5_value(pt_group, "HAS_SURFACE_DIFFUSION")
                config["has_surfDiff"] = "Yes" if has_surf_diff else "No"

    if resolution == "1D (radial coordinate)":
        _extract_particle_core_config(config, pt_group, prefix="")

    _extract_particle_geometry(config, pt_group, h5_unit_group, prefix="")

    return config


def get_reaction_type(group, phase):
    """Return the reaction model of the first <phase>_reaction_XXX subgroup, if any.

    This is the reaction interface of CADET-Core v6: a NREAC_<PHASE> count next to
    one subgroup per reaction, each naming its model in TYPE. The phase is one of
    "liquid", "solid" or "cross_phase".
    """

    if group is None:
        return None

    count = get_h5_value(group, "NREAC_" + phase.upper())
    if count is None or count < 1:
        return None

    reaction_group = group.get(f"{phase}_reaction_000")
    if reaction_group is None:
        return None

    return get_h5_value(reaction_group, "TYPE")


def _apply_reaction_model(config, cadet_reaction_model):
    """Store the generator's reaction model label, warning about unsupported models."""

    mapped = CADET_reaction_model_map.get(cadet_reaction_model, "Arbitrary")
    if mapped == "Arbitrary":
        st.sidebar.warning(
            f"Reaction model {cadet_reaction_model} not implemented in CADET-Equations, default to arbitrary reaction"
        )
    config["reaction_model"] = mapped


def _extract_reaction_config(config, h5_unit_group):
    """Extract the bulk liquid reaction configuration from an HDF5 unit group."""

    reaction_model_bulk = get_reaction_type(h5_unit_group, "liquid")

    if reaction_model_bulk is None:
        # pre-v6 interface: a single reaction model named at unit level
        reaction_bulk_group = h5_unit_group.get("reaction_bulk")
        if reaction_bulk_group is not None:
            reaction_model_bulk = get_h5_value(reaction_bulk_group, "REACTION_MODEL")

        if reaction_model_bulk is None:
            reaction_model_bulk = get_h5_value(h5_unit_group, "REACTION_MODEL")

    if reaction_model_bulk is not None and reaction_model_bulk != "NONE":
        config["has_reaction_bulk"] = "Yes"
        _apply_reaction_model(config, reaction_model_bulk)


def _extract_particle_reaction_config(config, pt_group):
    """Extract particle reaction configuration from a particle_type_XXX group.

    Cross-phase reactions are reported as solid phase reactions, since in CADET-Equations
    the particle solid reaction term is the one that depends on both phases.
    """

    liquid = get_reaction_type(pt_group, "liquid")
    solid = get_reaction_type(pt_group, "solid")
    cross_phase = get_reaction_type(pt_group, "cross_phase")

    if liquid is not None and liquid != "NONE":
        config["has_reaction_particle_liquid"] = "Yes"
        _apply_reaction_model(config, liquid)

    solid_phase = solid if solid is not None and solid != "NONE" else cross_phase
    if solid_phase is not None and solid_phase != "NONE":
        config["has_reaction_particle_solid"] = "Yes"
        _apply_reaction_model(config, solid_phase)

    if "has_reaction_particle_liquid" in config or "has_reaction_particle_solid" in config:
        # particle reactions are only offered in developer mode
        config["dev_mode"] = True
        config[_DEV_MODE_REQUIRED] = True


def _extract_particle_geometry(config, pt_group, h5_unit_group, prefix="particle_"):
    """Map PAR_GEOM to the generator's particle geometry.

    PAR_GEOM moved from the unit's discretization group into particle_type_XXX.
    """

    par_geom = get_h5_value(pt_group, "PAR_GEOM") if pt_group is not None else None

    if par_geom is None:
        disc_group = h5_unit_group.get("discretization") if h5_unit_group is not None else None
        par_geom = get_h5_value(disc_group, "PAR_GEOM") if disc_group is not None else None

    if par_geom is None:
        return

    geometry = CADET_particle_geometry_map.get(par_geom)
    if geometry is None:
        st.sidebar.warning(f"Particle geometry {par_geom} is not available in CADET-Equations, assuming a sphere.")
        return

    if geometry != "Sphere":
        # non-spherical particles are only offered in developer mode
        config[prefix + "geometry"] = geometry
        config["dev_mode"] = True
        config[_DEV_MODE_REQUIRED] = True


def _extract_particle_core_config(config, group, prefix="particle_"):
    """Extract particle core radius config. Shared between v5 (unit group) and v6 (particle_type group)."""
    config[prefix + "has_core"] = "No"
    parCore = get_h5_value(group, "PAR_CORERADIUS")
    if parCore is not None:
        if parCore > 0.0:
            config[prefix + "has_core"] = "Yes"
            config["advanced_mode"] = "On"


CADET_crystallization_unit_types = ["CSTR", "LUMPED_RATE_MODEL_WITHOUT_PORES"]

AGGREGATION_KERNEL_MAP = {
    0: "Constant",
    1: "Brownian",
    2: "Smoluchowski",
    3: "Golovin",
    4: "Differential force",
}


def _is_crystallization_unit(h5_unit_group):
    if get_reaction_type(h5_unit_group, "liquid") == "CRYSTALLIZATION":
        return True
    reaction_model = get_h5_value(h5_unit_group, "REACTION_MODEL")
    return reaction_model == "CRYSTALLIZATION"


def _get_crystallization_reaction_group(h5_unit_group):
    for name in ("liquid_reaction_000", "reaction_bulk", "reaction"):
        grp = h5_unit_group.get(name)
        if grp is not None:
            return grp
    return None


def extract_crystallization_config(unit_type, h5_unit_group):

    config = {}
    config["model_type"] = "Crystallization"
    config["dev_mode"] = True
    config["var_format"] = "CADET"
    config["show_eq_description"] = True
    config["model_assumptions"] = True

    if unit_type == "CSTR":
        config["cry_column_type"] = "CSTR"
    else:
        config["cry_column_type"] = "DPFR"

    reaction_group = _get_crystallization_reaction_group(h5_unit_group)

    if reaction_group is None:
        config["cry_has_primary_formation"] = "Yes"
        config["cry_has_aggregation"] = "No"
        config["cry_has_fragmentation"] = "No"
        return config

    cry_mode = get_h5_value(reaction_group, "CRY_MODE")
    if cry_mode is None:
        cry_mode = 1

    has_primary = bool(cry_mode & 1)
    has_agg = bool(cry_mode & 2)
    has_frag = bool(cry_mode & 4)

    config["cry_has_primary_formation"] = "Yes" if has_primary else "No"

    if has_primary:
        gd_rate = get_h5_value(reaction_group, "CRY_GROWTH_DISPERSION_RATE")
        config["cry_has_growth_dispersion"] = "Yes" if (gd_rate is not None and gd_rate > 0) else "No"

        sec_rate = get_h5_value(reaction_group, "CRY_SECONDARY_NUCLEATION_RATE")
        config["cry_has_secondary_nucleation"] = "Yes" if (sec_rate is not None and sec_rate > 0) else "No"

        cry_p = get_h5_value(reaction_group, "CRY_P")
        config["cry_size_dependent_growth"] = "Yes" if (cry_p is not None and cry_p != 0) else "No"

    config["cry_has_aggregation"] = "Yes" if has_agg else "No"

    if has_agg:
        agg_idx = get_h5_value(reaction_group, "CRY_AGGREGATION_INDEX")
        if agg_idx is None:
            agg_idx = 0
        config["cry_aggregation_kernel"] = AGGREGATION_KERNEL_MAP.get(agg_idx, "Constant")

    config["cry_has_fragmentation"] = "Yes" if has_frag else "No"

    if config["cry_column_type"] == "DPFR":
        col_disp = get_h5_value(h5_unit_group, "COL_DISPERSION")
        config["cry_has_axial_dispersion"] = "Yes" if (col_disp is not None and col_disp > 1e-20) else "No"

    return config


def get_config_from_CADET_h5(h5_filename, unit_idx):

    with h5py.File(h5_filename, "r") as f:
        model_group = f["input/model"]

        if unit_idx == "-01":
            unit_keys = [k for k in model_group if re.match(r"^unit_\d{3}$", k)]
            unit_keys.sort(key=lambda x: int(x.split("_")[1]))

            for unit_key in unit_keys:
                unit_group = model_group[unit_key]

                unit_type = get_h5_value(unit_group, "UNIT_TYPE")

                if unit_type is not None:
                    if unit_type in CADET_column_unit_types:
                        st.sidebar.success(
                            unit_type + " was found in " + re.sub(r"input/model/", "", unit_key) + " and is applied!"
                        )

                        if _is_crystallization_unit(unit_group):
                            return extract_crystallization_config(unit_type, unit_group)

                        return extract_config_data_from_unit(unit_type, unit_group)

            st.sidebar.error("No supported column unit type was found in the file.")

        else:
            if "unit_" + unit_idx in model_group:
                unit_group = model_group["unit_" + unit_idx]
                unit_type = get_h5_value(unit_group, "UNIT_TYPE")
            else:
                st.sidebar.error(f"unit_{unit_idx} does not exist!")
                return None

            if unit_type in CADET_column_unit_types:
                st.sidebar.success(f"Unit type {unit_type} is applied.")

                if _is_crystallization_unit(unit_group):
                    return extract_crystallization_config(unit_type, unit_group)

                return extract_config_data_from_unit(unit_type, unit_group)

            else:
                st.sidebar.error(f"Equations for {unit_type} are not available.")

        return None
