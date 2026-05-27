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


@dataclass
class Crystallization:

    var_format: Literal
    column_type: str = "CSTR"
    cry_mode: int = 1
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

        mode_labels = list(eq.CRY_MODES.values())
        mode_keys = list(eq.CRY_MODES.keys())
        selected_label = st.sidebar.selectbox(
            "Crystallization mode (CRY_MODE)",
            mode_labels, key=r"cry_mode")
        self.cry_mode = mode_keys[mode_labels.index(selected_label)]

        if eq.cry_has_primary_formation(self.cry_mode):
            with st.sidebar.expander("Growth and nucleation", expanded=True):
                self.size_dependent_growth = st.selectbox(
                    "Size-dependent growth", ["No", "Yes"],
                    key=r"cry_size_dependent_growth") == "Yes"
                self.has_growth_dispersion = st.selectbox(
                    "Add growth dispersion", ["No", "Yes"],
                    key=r"cry_has_growth_dispersion") == "Yes"
                self.has_secondary_nucleation = st.selectbox(
                    "Add secondary nucleation", ["No", "Yes"],
                    key=r"cry_has_secondary_nucleation") == "Yes"

        if eq.cry_has_aggregation(self.cry_mode):
            with st.sidebar.expander("Aggregation", expanded=True):
                kernel_labels = list(eq.AGGREGATION_KERNELS.values())
                kernel_keys = list(eq.AGGREGATION_KERNELS.keys())
                selected_kernel = st.selectbox(
                    "Aggregation kernel", kernel_labels,
                    key=r"cry_aggregation_kernel")
                self.aggregation_kernel_index = kernel_keys[kernel_labels.index(selected_kernel)]

        self.fill_vars_and_params()

    def model_name(self):
        name = "Crystallization"
        if self.column_type == "CSTR":
            name += " in CSTR"
        else:
            name += " in DPFR"

        parts = []
        if eq.cry_has_primary_formation(self.cry_mode):
            parts.append("growth/nucleation")
        if eq.cry_has_aggregation(self.cry_mode):
            parts.append("aggregation")
        if eq.cry_has_fragmentation(self.cry_mode):
            parts.append("fragmentation")

        if parts:
            name += " with " + ", ".join(parts)
        return name

    def model_assumptions(self):
        return {
            "Model assumptions": eq.cry_assumptions(self.column_type, self.cry_mode)
        }

    def available_CADET_Core(self):
        return 1

    def available_CADET_Process(self):
        return -1

    def fill_vars_and_params(self):
        has_primary = eq.cry_has_primary_formation(self.cry_mode)
        has_agg = eq.cry_has_aggregation(self.cry_mode)
        has_frag = eq.cry_has_fragmentation(self.cry_mode)

        vp = [
            {"Group": 0, "Symbol": r"t", "Description": "time coordinate",
             "Unit": r"s", "Dependence": r"\text{independent variable}"},
            {"Group": 0, "Symbol": r"x", "Description": "particle size (internal coordinate)",
             "Unit": r"m", "Dependence": r"\text{independent variable}"},
            {"Group": 1, "Symbol": r"n", "Description": "particle number density",
             "Unit": r"\frac{1}{m \cdot m^3}", "Dependence": r"t, x"},
            {"Group": 1, "Symbol": r"c", "Description": "solute concentration",
             "Unit": r"\frac{kg}{m^3}", "Dependence": r"t"},
        ]

        if self.column_type == "CSTR":
            vp.append({"Group": 1.5, "Symbol": r"V", "Description": "reactor volume",
                        "Unit": r"m^3", "Dependence": r"t"})
            vp.append({"Group": 2, "Symbol": r"F_{\mathrm{in}}", "Description": "volumetric inflow rate",
                        "Unit": r"\frac{m^3}{s}", "Dependence": r"\text{constant}"})
            vp.append({"Group": 2, "Symbol": r"F_{\mathrm{out}}", "Description": "volumetric outflow rate",
                        "Unit": r"\frac{m^3}{s}", "Dependence": r"\text{constant}"})
        else:
            vp.append({"Group": 0, "Symbol": r"z", "Description": "axial coordinate",
                        "Unit": r"m", "Dependence": r"\text{independent variable}"})
            vp[2]["Dependence"] = r"t, x, z"
            vp[3]["Dependence"] = r"t, z"
            vp.append({"Group": 2, "Symbol": r"v_{\mathrm{ax}}", "Description": "axial velocity",
                        "Unit": r"\frac{m}{s}", "Dependence": r"\text{constant}"})
            vp.append({"Group": -1, "Symbol": r"L", "Description": "reactor length",
                        "Unit": r"m", "Dependence": r"\text{constant}"})
            if self.has_axial_dispersion:
                vp.append({"Group": 6, "Symbol": r"D_{\mathrm{ax}}", "Description": "axial dispersion coefficient",
                            "Unit": r"\frac{m^2}{s}", "Dependence": r"\text{constant}"})

        if has_primary:
            vp.append({"Group": 3, "Symbol": r"v_G", "Description": "growth rate",
                        "Unit": r"\frac{m}{s}", "Dependence": r"s, x"})
            vp.append({"Group": 3, "Symbol": r"s", "Description": "relative supersaturation",
                        "Unit": r"-", "Dependence": r"c"})
            vp.append({"Group": 3, "Symbol": r"c_{\mathrm{eq}}", "Description": "equilibrium solubility",
                        "Unit": r"\frac{kg}{m^3}", "Dependence": r"\text{constant}"})
            vp.append({"Group": 3, "Symbol": r"B_0", "Description": "total nucleation rate",
                        "Unit": r"\frac{1}{m^3 \cdot s}", "Dependence": r"s"})
            vp.append({"Group": 4, "Symbol": r"k_g", "Description": "growth rate constant",
                        "Unit": r"\frac{m}{s}", "Dependence": r"\text{constant}"})
            vp.append({"Group": 4, "Symbol": r"g", "Description": "growth exponent",
                        "Unit": r"-", "Dependence": r"\text{constant}"})
            vp.append({"Group": 4, "Symbol": r"a", "Description": "growth parameter",
                        "Unit": r"-", "Dependence": r"\text{constant}"})
            if self.size_dependent_growth:
                vp.append({"Group": 4, "Symbol": r"\gamma", "Description": "size-dependence quantifier",
                            "Unit": r"-", "Dependence": r"\text{constant}"})
                vp.append({"Group": 4, "Symbol": r"p", "Description": "size-dependence exponent",
                            "Unit": r"-", "Dependence": r"\text{constant}"})
            if self.has_growth_dispersion:
                vp.append({"Group": 6, "Symbol": r"D_g", "Description": "growth dispersion rate",
                            "Unit": r"\frac{m^2}{s}", "Dependence": r"\text{constant}"})
            vp.append({"Group": 5, "Symbol": r"k_p", "Description": "primary nucleation rate constant",
                        "Unit": r"\frac{1}{m^3 \cdot s}", "Dependence": r"\text{constant}"})
            vp.append({"Group": 5, "Symbol": r"u", "Description": "primary nucleation exponent",
                        "Unit": r"-", "Dependence": r"\text{constant}"})
            if self.has_secondary_nucleation:
                vp.append({"Group": 5, "Symbol": r"k_b", "Description": "secondary nucleation rate constant",
                            "Unit": r"\frac{1}{m^3 \cdot s}", "Dependence": r"\text{constant}"})
                vp.append({"Group": 5, "Symbol": r"b", "Description": "secondary nucleation exponent",
                            "Unit": r"-", "Dependence": r"\text{constant}"})
                vp.append({"Group": 5, "Symbol": r"k", "Description": "suspension density exponent",
                            "Unit": r"-", "Dependence": r"\text{constant}"})
                vp.append({"Group": 5, "Symbol": r"M", "Description": "suspension density",
                            "Unit": r"\frac{kg}{m^3}", "Dependence": r"n"})
            vp.append({"Group": 5, "Symbol": r"x_c", "Description": "critical (minimum) particle size",
                        "Unit": r"m", "Dependence": r"\text{constant}"})
            vp.append({"Group": -1, "Symbol": r"\rho", "Description": "nuclei mass density",
                        "Unit": r"\frac{kg}{m^3}", "Dependence": r"\text{constant}"})
            vp.append({"Group": -1, "Symbol": r"k_v", "Description": "volumetric shape factor",
                        "Unit": r"-", "Dependence": r"\text{constant}"})

        if has_agg:
            vp.append({"Group": 7, "Symbol": r"\beta", "Description": "aggregation kernel",
                        "Unit": r"\frac{m^3}{s}", "Dependence": r"x, \lambda"})
            vp.append({"Group": 7, "Symbol": r"\beta_0", "Description": "aggregation rate constant",
                        "Unit": r"\frac{m^3}{s}", "Dependence": r"\text{constant}"})

        if has_frag:
            vp.append({"Group": 8, "Symbol": r"S(x)", "Description": "selection (fragmentation rate) function",
                        "Unit": r"\frac{1}{s}", "Dependence": r"x"})
            vp.append({"Group": 8, "Symbol": r"b(x \mid \lambda)", "Description": "breakage probability density",
                        "Unit": r"\frac{1}{m}", "Dependence": r"x, \lambda"})
            vp.append({"Group": 8, "Symbol": r"S_0", "Description": "fragmentation rate constant",
                        "Unit": r"\frac{1}{s}", "Dependence": r"\text{constant}"})
            vp.append({"Group": 8, "Symbol": r"\alpha", "Description": "fragmentation exponent",
                        "Unit": r"-", "Dependence": r"\text{constant}"})
            vp.append({"Group": 8.1, "Symbol": r"\gamma_f", "Description": "daughter distribution parameter",
                        "Unit": r"-", "Dependence": r"\text{constant}"})

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
