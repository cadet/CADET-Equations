# -*- coding: utf-8 -*-
"""
This script implements auxiliary models
"""

from attr import dataclass
from typing import List
import re

from src.utils import format_variables


@dataclass
class AuxiliaryModel:

    aux_model: str
    aux_model_list = ['INLET', 'OUTLET']
    N_c: int = -1
    N_p:int = 0
    
    var_format: str = "CADET"
    unit_system: str = "SI"

    def __post_init__(self):
        
        if self.aux_model not in self.aux_model_list:
            raise ValueError(
                f"Invalid auxiliary model: {self.aux_model}. Must be one of {self.aux_model_list}.")
        
        self.vars_and_params()

    def vars_and_params(self):
        
        if self.aux_model == 'INLET':
            state_deps = r"t; i"
            param_deps = r"t; i"
    
            self.vars_and_params = [
                {"Group" : 0, "Symbol": r"t", "Description": r"time coordinate", "Unit": r"s", "Dependence": r"\text{independent variable}", "Property": r"\in (0, T^{\mathrm{end}})"},
                {"Group" : 1, "Symbol": r"c^{\b, in}_i", "Description": r"feed of bulk liquid concentration", "Unit": r"\frac{mol}{m^3}", "Dependence" : state_deps, "Domain": "(0, T^{\mathrm{end}})"},
                {"Group" : -1, "Symbol": r"T^{\mathrm{end}}", "Description": r"process end time", "Unit": r"s", "Dependence": r"\text{constant}", "Property": r" > 0"},
                {"Group" : -0.99, "Symbol": r"0 \leq t_1 < \dots < t_{N^{\mathrm{sec}}+1} \leq T^{\mathrm{end}}", "Description": r"decomposition of the simulation time interval", "Unit": r"s", "Dependence": r"\text{constant}"},
                {"Group" : -0.2, "Symbol": r"N^\mathrm{sec}", "Description": r"number of time sections", "Unit": r"-", "Dependence": r"-", "Property": r"\in \mathbb{N}"},
                {"Group" : -0.2, "Symbol": r"N^\mathrm{c}", "Description": r"number of components", "Unit": r"-", "Dependence": r"-", "Property": r"\in \mathbb{N}"},
                {"Group" : -0.1, "Symbol": r"i", "Description": r"component index", "Unit": r"s", "Dependence": r"-", "Property": r"-"},
                ]
            
            for var_ in self.vars_and_params:
                var_["Symbol"] = format_variables(var_["Symbol"], self.var_format)
                
            self.vars_and_params = sorted(self.vars_and_params, key=lambda x: x['Group'])
        else:
            self.vars_and_params = []
        
    def model_name(self):

        return self.aux_model

    def model_assumptions(self):
        
        if self.aux_model == 'INLET':
            return {
                "General model assumptions": [
                    r"An inlet unit operation defines a source of inflow (feed) to the system.",
					r"The feed of each component can be described by a piecewise cubic polynomial."
				]
            }
        elif self.aux_model == 'OUTLET':
            return {
                "General model assumptions": [
					r"An outlet unit operation represents a sink in the system, describing outflow from the system."
                ]
            }
        
        return None

    def model_equation(self):
        
        if self.aux_model == 'INLET':
            return r"""
\begin{equation}
    c_i^{\b, in} = \sum_{k=1}^{N^{\mathrm{sec}}} \mathbb{R}_{[t_k, t_{k+1})}(t) \left[a_{k,i} (t-t_k)^3 + b_{k,i} (t-t_k)^2 + d_{k,i} (t-t_k) + f_{k,i} \right].
\end{equation}
"""
        else:
            return None

    def vars_params_description(self):

        description_ = ""

        idx_ = 1
        num_VP = len(self.vars_and_params)

        for thing in self.vars_and_params:

            if thing.get("Group", -1) < 0: # dont print symbols with negative group no.
                num_VP -= 1
                continue

            if not idx_ == 1:
                description_ += ", " if idx_ < num_VP else ", and "
            description_ += r"$" + thing["Symbol"]

            if not thing.get("Domain", "-") == "-":
                description_ += r"\colon " + re.sub(r"\$", "", thing["Domain"]) + r" \mapsto \mathbb{R}"

            description_ += thing.get("Property", "") + r"$"
            
            description_ += " is the " + thing["Description"]

            idx_ += 1
                
        return description_ + "."