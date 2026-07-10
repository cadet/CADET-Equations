# -*- coding: utf-8 -*-
"""
This script implements the `Crystallization` dataclass for population balance
models used in crystallization / precipitation processes.
"""

from dataclasses import dataclass, field
from typing import List, Literal, Optional
import streamlit as st

from src import equations as eq
from src.utils import format_variables
from src.units import get_unit


@dataclass
class Crystallization:

    var_format: Literal
    unit_system: str = "SI"
    column_type: str = "CSTR"
    has_primary_formation: bool = True
    has_aggregation: bool = False
    has_fragmentation: bool = False
    has_axial_dispersion: bool = False
    has_growth_dispersion: bool = False
    size_dependent_growth: bool = False
    has_secondary_nucleation: bool = False
    aggregation_kernel_index: int = 0

    vars_and_params: List[dict] = field(default_factory=list)

    def __post_init__(self):

        st.sidebar.write("Configure crystallization model")

        self.column_type = st.sidebar.selectbox(
            "Reactor type", ["CSTR", "DPFR"],
            key=r"cry_column_type")

        if self.column_type == "DPFR":
            self.has_axial_dispersion = st.sidebar.selectbox(
                "Add axial dispersion", ["Yes", "No"],
                key=r"cry_has_axial_dispersion") == "Yes"

        st.sidebar.write("Crystallization mechanisms")
        self.has_primary_formation = st.sidebar.selectbox(
            "Primary particle formation (nucleation and growth)",
            ["Yes", "No"], key=r"cry_has_primary_formation") == "Yes"

        if self.has_primary_formation:
            with st.sidebar.expander("Nucleation and growth", expanded=True):
                self.size_dependent_growth = st.selectbox(
                    "Size-dependent growth", ["No", "Yes"],
                    key=r"cry_size_dependent_growth") == "Yes"
                self.has_growth_dispersion = st.selectbox(
                    "Add growth dispersion", ["No", "Yes"],
                    key=r"cry_has_growth_dispersion") == "Yes"
                self.has_secondary_nucleation = st.selectbox(
                    "Add secondary nucleation", ["No", "Yes"],
                    key=r"cry_has_secondary_nucleation") == "Yes"

        self.has_aggregation = st.sidebar.selectbox(
            "Aggregation",
            ["No", "Yes"], key=r"cry_has_aggregation") == "Yes"

        if self.has_aggregation:
            with st.sidebar.expander("Aggregation", expanded=True):
                kernel_labels = list(eq.AGGREGATION_KERNELS.values())
                kernel_keys = list(eq.AGGREGATION_KERNELS.keys())
                selected_kernel = st.selectbox(
                    "Aggregation kernel", kernel_labels,
                    key=r"cry_aggregation_kernel")
                self.aggregation_kernel_index = kernel_keys[kernel_labels.index(selected_kernel)]

        self.has_fragmentation = st.sidebar.selectbox(
            "Fragmentation",
            ["No", "Yes"], key=r"cry_has_fragmentation") == "Yes"

        if not self.has_primary_formation and not self.has_aggregation and not self.has_fragmentation:
            st.sidebar.warning("At least one mechanism must be enabled.")
            self.has_primary_formation = True

        self.fill_vars_and_params()

    def model_name(self):
        name = "Crystallization"
        if self.column_type == "CSTR":
            name += " in CSTR"
        else:
            name += " in DPFR"

        parts = []
        if self.has_primary_formation:
            parts.append("nucleation and growth")
        if self.has_aggregation:
            parts.append("aggregation")
        if self.has_fragmentation:
            parts.append("fragmentation")

        if parts:
            name += " with " + ", ".join(parts)
        return name

    def aggregation_kernel_name(self):
        """Return the human-readable aggregation kernel label."""
        return eq.AGGREGATION_KERNELS.get(self.aggregation_kernel_index, eq.AGGREGATION_KERNELS[0])

    def model_assumptions(self):
        return {
            "Model assumptions": eq.cry_assumptions(
                self.column_type, self.has_primary_formation,
                self.has_aggregation, self.has_fragmentation)
        }

    def available_CADET_Core(self):
        return 1

    def available_CADET_Process(self):
        return -1

    def available_CADET_SemiAnalytic(self):
        return -1

    def fill_vars_and_params(self):
        has_primary = self.has_primary_formation
        has_agg = self.has_aggregation
        has_frag = self.has_fragmentation
        u = lambda key: get_unit(key, self.unit_system)

        vp = [
            {"Group": 0, "Symbol": r"t", "Description": "time coordinate",
             "Unit": u("time"), "Dependence": r"\text{independent variable}"},
            {"Group": 0, "Symbol": r"x", "Description": "particle size (internal coordinate)",
             "Unit": u("length"), "Dependence": r"\text{independent variable}"},
            {"Group": 1, "Symbol": r"n", "Description": "particle number density",
             "Unit": u("number_density"), "Dependence": r"t, x"},
            {"Group": 1, "Symbol": r"c", "Description": "solute concentration",
             "Unit": u("concentration_mass"), "Dependence": r"t"},
        ]

        if self.column_type == "CSTR":
            vp.append({"Group": 1.5, "Symbol": r"V", "Description": "reactor volume",
                        "Unit": u("volume"), "Dependence": r"t"})
            vp.append({"Group": 2, "Symbol": r"F_{\mathrm{in}}", "Description": "volumetric inflow rate",
                        "Unit": u("volumetric_flow"), "Dependence": r"\text{constant}"})
            vp.append({"Group": 2, "Symbol": r"F_{\mathrm{out}}", "Description": "volumetric outflow rate",
                        "Unit": u("volumetric_flow"), "Dependence": r"\text{constant}"})
        else:
            vp.append({"Group": 0, "Symbol": r"z", "Description": "axial coordinate",
                        "Unit": u("length"), "Dependence": r"\text{independent variable}"})
            vp[2]["Dependence"] = r"t, x, z"
            vp[3]["Dependence"] = r"t, z"
            vp.append({"Group": 2, "Symbol": r"v_{\mathrm{ax}}", "Description": "axial velocity",
                        "Unit": u("velocity"), "Dependence": r"\text{constant}"})
            vp.append({"Group": -1, "Symbol": r"L", "Description": "reactor length",
                        "Unit": u("length"), "Dependence": r"\text{constant}"})
            if self.has_axial_dispersion:
                vp.append({"Group": 6, "Symbol": r"D_{\mathrm{ax}}", "Description": "axial dispersion coefficient",
                            "Unit": u("diffusion"), "Dependence": r"\text{constant}"})

        if has_primary:
            vp.append({"Group": 3, "Symbol": r"v_G", "Description": "growth rate",
                        "Unit": u("velocity"), "Dependence": r"s, x"})
            vp.append({"Group": 3, "Symbol": r"s", "Description": "relative supersaturation",
                        "Unit": u("dimensionless"), "Dependence": r"c"})
            vp.append({"Group": 3, "Symbol": r"c_{\mathrm{eq}}", "Description": "equilibrium solubility",
                        "Unit": u("concentration_mass"), "Dependence": r"\text{constant}"})
            vp.append({"Group": 3, "Symbol": r"B_0", "Description": "total nucleation rate",
                        "Unit": u("nucleation_rate"), "Dependence": r"s"})
            vp.append({"Group": 4, "Symbol": r"k_g", "Description": "growth rate constant",
                        "Unit": u("velocity"), "Dependence": r"\text{constant}"})
            vp.append({"Group": 4, "Symbol": r"g", "Description": "growth exponent",
                        "Unit": u("dimensionless"), "Dependence": r"\text{constant}"})
            vp.append({"Group": 4, "Symbol": r"a", "Description": "growth parameter",
                        "Unit": u("dimensionless"), "Dependence": r"\text{constant}"})
            if self.size_dependent_growth:
                vp.append({"Group": 4, "Symbol": r"\gamma", "Description": "size-dependence quantifier",
                            "Unit": u("dimensionless"), "Dependence": r"\text{constant}"})
                vp.append({"Group": 4, "Symbol": r"p", "Description": "size-dependence exponent",
                            "Unit": u("dimensionless"), "Dependence": r"\text{constant}"})
            if self.has_growth_dispersion:
                vp.append({"Group": 6, "Symbol": r"D_g", "Description": "growth dispersion rate",
                            "Unit": u("diffusion"), "Dependence": r"\text{constant}"})
            vp.append({"Group": 5, "Symbol": r"k_p", "Description": "primary nucleation rate constant",
                        "Unit": u("nucleation_rate"), "Dependence": r"\text{constant}"})
            vp.append({"Group": 5, "Symbol": r"u", "Description": "primary nucleation exponent",
                        "Unit": u("dimensionless"), "Dependence": r"\text{constant}"})
            if self.has_secondary_nucleation:
                vp.append({"Group": 5, "Symbol": r"k_b", "Description": "secondary nucleation rate constant",
                            "Unit": u("nucleation_rate"), "Dependence": r"\text{constant}"})
                vp.append({"Group": 5, "Symbol": r"b", "Description": "secondary nucleation exponent",
                            "Unit": u("dimensionless"), "Dependence": r"\text{constant}"})
                vp.append({"Group": 5, "Symbol": r"k", "Description": "suspension density exponent",
                            "Unit": u("dimensionless"), "Dependence": r"\text{constant}"})
                vp.append({"Group": 5, "Symbol": r"M", "Description": "suspension density",
                            "Unit": u("concentration_mass"), "Dependence": r"n"})
            vp.append({"Group": 5, "Symbol": r"x_c", "Description": "critical (minimum) particle size",
                        "Unit": u("length"), "Dependence": r"\text{constant}"})
            vp.append({"Group": -1, "Symbol": r"\rho", "Description": "nuclei mass density",
                        "Unit": u("concentration_mass"), "Dependence": r"\text{constant}"})
            vp.append({"Group": -1, "Symbol": r"k_v", "Description": "volumetric shape factor",
                        "Unit": u("dimensionless"), "Dependence": r"\text{constant}"})

        if has_agg:
            vp.append({"Group": 7, "Symbol": r"\beta", "Description": "aggregation kernel",
                        "Unit": u("volumetric_flow"), "Dependence": r"x, \lambda"})
            vp.append({"Group": 7, "Symbol": r"\beta_0", "Description": "aggregation rate constant",
                        "Unit": u("volumetric_flow"), "Dependence": r"\text{constant}"})

        if has_frag:
            vp.append({"Group": 8, "Symbol": r"S(x)", "Description": "selection (fragmentation rate) function",
                        "Unit": u("rate_first_order"), "Dependence": r"x"})
            vp.append({"Group": 8, "Symbol": r"b(x \mid \lambda)", "Description": "breakage probability density",
                        "Unit": u("inverse_length"), "Dependence": r"x, \lambda"})
            vp.append({"Group": 8, "Symbol": r"S_0", "Description": "fragmentation rate constant",
                        "Unit": u("rate_first_order"), "Dependence": r"\text{constant}"})
            vp.append({"Group": 8, "Symbol": r"\alpha", "Description": "fragmentation exponent",
                        "Unit": u("dimensionless"), "Dependence": r"\text{constant}"})
            vp.append({"Group": 8.1, "Symbol": r"\gamma_f", "Description": "daughter distribution parameter",
                        "Unit": u("dimensionless"), "Dependence": r"\text{constant}"})

        for var_ in vp:
            var_["Symbol"] = format_variables(var_["Symbol"], self.var_format)

        self.vars_and_params = sorted(vp, key=lambda v: v["Group"])

    def vars_params_description(self):
        description_ = ""
        idx_ = 1
        num_VP = len(self.vars_and_params)
        for thing in self.vars_and_params:
            if thing.get("Group", -1) < 0:
                num_VP -= 1
                continue
            if idx_ != 1:
                description_ += ", " if idx_ < num_VP else ", and "
            description_ += r"$" + thing["Symbol"]
            description_ += thing.get("Property", "") + r"$"
            description_ += " is the " + thing["Description"]
            idx_ += 1
        return description_ + "."
