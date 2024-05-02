# -*- coding: utf-8 -*-

"""
This is the interface for the TREND 4.1 dynamic library call.
Expose the so/dll functionality for two supported systems {LINUX}.
"""
__author__ = 'David Celný'
__version__ = '0.0.1'

import sys
import os
import re
import subprocess
from contextlib import contextmanager
from typing import Optional
# own files
from eos_abc import Eos
from TREND_4.helper import Property, Substances
# from setup import *
# import calculation

# @contextmanager
# def suppress_stdout():
#     with open(os.devnull, "w") as devnull:
#         old_stdout = sys.stdout
#         sys.stdout = devnull
#         try:  
#             yield
#         finally:
#             sys.stdout = old_stdout

def get_result(capture:str)-> float:
    keyword = 'result:'
    cut = capture.find(keyword)
    cut2 = capture.find('errorflag:')
    # print(capture[cut+len(keyword):cut2]) #DEBUG
    if cut > -1:
        return float(capture[cut+len(keyword):cut2])
    return None

# trend interface related functionality
class trend_eos:
    name_tag = "trend"
    programpath = os.path.realpath("change_this_path_to_point_at_trend_main_program_executable")
    def _setup_interface_forms(self, eos_type):
        self.eos_indicator = f"           {eos_type}           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0\n"     
        self.mix_indicator = "           1\n"
        self.capture_status = -1

    def __init__(self, eos_type=6, lib_path_str: Optional[str] = None):
        self.name_tag += f"_{eos_type}"
        self.inputpath = os.path.abspath("input.txt")
        self._setup_interface_forms(eos_type)
        # self.calculation = None
        return
    def __str__(self):
        tmp_str =  f"EOS        : {self.eostype}\n"
        tmp_str += f"path       : {self.dll_path}\n"
        methods = [method for method in dir(trend_eos) if callable(getattr(trend_eos, method)) and not method.startswith("_")]
        tmp_str += f"methods    : {methods},\n"
        return tmp_str

    def __eq__(self, other):
        if (self.programpath != other.programpath):
            return False
        # if (self.inputpath != other.inputpath):
        #   return False
        if (self.eos_indicator != other.eos_indicator):
            return False
        # if (self.mix_indicator != other.mix_indicator):
        #   return False
        if (self.capture_status != other.capture_status):
            return False
        return True

    def __ne__(self, other):
        return not(self == other)

    def _call(self, calc_target, input_target, prop1, prop2, substance_names, substance_molrat, is_property_call:bool=False):
        # facilitate the call to the dll and return the unparsed output and status
        # mode has to be stringized number "0" for EoS and "1" for property
        substances = Substances(substance_names, substance_molrat)
        with open(self.inputpath,mode="w") as file:
            file.write(f"{calc_target}\n")
            file.write(f"{input_target}\n")
            file.write(f"{prop1:f}\n")
            file.write(f"{prop2:f}\n")
            file.write(substances.get_trend_name_str())
            file.write(substances.get_trend_ratio_str())
            file.write(self.eos_indicator)
            file.write(self.mix_indicator)
            file.write(f"{int(is_property_call)}")

        # print([str(self.programpath), str(self.inputpath)])
        capture = subprocess.run(args=[str(self.programpath), str(self.inputpath)],
                                 capture_output=True,
                                 text=True)
        # print(capture.stdout) # print all the output - loads of stuff when debugging
        # TODO should detect errors here
        # capture_status = int(mode)
        return capture.stdout

    # def call_n(self, calc_target, input_target, prop1, prop2, substance_names, substance_molrat, verbose=False):
    #     capture = self._call(calc_target, input_target, prop1, prop2, substance_names, substance_molrat, False)

    #     # capture the last A_res output to ommit the long output
    #     # print(capture)
    #     last_section_capture = capture[capture.rfind("=== ===")+len("=== ==="):]
        
    #     if verbose:
    #         print(last_section_capture)
    #     res = re.findall(r"result:\s* ([-\d\.]*)", last_section_capture)
    #     # print(res)
    #     if len(res) == 0:
    #         print("Error: no result found in captured output")
    #         return None
    #     else:
    #         return float(res[0])


    def get_p(self, temperature, density, composition, substance, is_silent: bool= True) ->float:
        """ return the pressure in Pa """
        result = []
        for calc_type in ("p","DPDD","D2PDD2"):
        # for calc_type in ("p"):
            # density *= 1e3
            capture = self._call(calc_type, "TD", temperature, density, substance, composition, is_property_call=False)
            # print("capture in get_p: \n",capture) #DEBUG
            # capture the last A_res output to ommit the long output
            last_section_capture = capture[capture.rfind("=== ===")+len("=== ==="):]
            
            if not is_silent:
                print(last_section_capture)
                print(type(last_section_capture))

            res = get_result(last_section_capture)
            # print(res) #DEBUG

            if res is None:
                print(f"Error: no result found in captured output for {calc_type}")
            result.append(res) 

        #NOTE tested to be working with these conversion rules 07.05.2022 vs nist for methane
        result[0] *= 1e6 #NOTE conversion to units Pa
        result[1] *= 1e6 #NOTE conversion to units Pa*m**3/mol
        result[2] *= 1e6 #NOTE conversion to units Pa*m**3/mol
        return tuple(result)

    def get_mu_fug(self, temperature, density, composition, substance, is_silent: bool= True) ->float:
        """ return the residual chemical potential for pure fluid only """
        capture = self._call("CPOTR", "TD", temperature, density, substance, composition, is_property_call=False)

        # capture the last A_res output to ommit the long output
        last_section_capture = capture[capture.rfind("=== ===")+len("=== ==="):]
        
        if not is_silent:
            print(last_section_capture)

        res = re.findall(r"result:\s* ([-\d\.]*)", last_section_capture)
        # print(res) #DEBUG

        if len(res) == 0:
            print("Error: no result found in captured output")
            return None, None
        else:
            return float(res[0]), None

    # def call_helmholtz(self, calculation, verbose=False):
    #     self.eval(calculation, "0")
        
    #     # capture the last A_res output to ommit the long output
    #     last_section_capture = self.capture[self.capture.rfind("=== ===")+len("=== ==="):]
        
    #     # print(last_section_capture) # print ontly the last of A_res
    #     helmoltz_res = re.findall(r"(A\w*)\s*: ((?:\s{1,8}[-E\d.]*){15})", last_section_capture)        
    #     if verbose:
    #         print("helmholtz results")
    #     for el_id,el in enumerate(helmoltz_res):
    #         tmp = [float(num) for num in el[1].split()]
    #         helmoltz_res[el_id] = (el[0].lower(),tmp)
    #         if verbose:
    #             print(helmoltz_res[el_id])
    #     if verbose:
    #         print(" ")

    #     return dict(helmoltz_res)

    # def call_special(self, calculation, verbose=False):
    #     self.eval(calculation, "0")

    #     # print(self.capture)
    #     # handling of the specially requested outputs
    #     special_res = re.findall(r"===== (\w*)\s*:\s*([-\d\.E]*)\s*=====", self.capture)
    #     # print(special_res)
    #     if verbose:
    #         print("requested inputs")
    #     for el_id,el in enumerate(special_res):
    #         tmp = [float(num) for num in el[1].split()]
    #         special_res[el_id] = (el[0],tmp)
    #         if verbose:
    #             print(special_res[el_id])
    #     if verbose:
    #         print(" ")

    #     return special_res

    # def call_prop(self, calculation, verbose=False):
    #     # used for obtaining only the properties in following order:
    #     # Mol weight, T_tripple, P-tripple, T_crit, P_crit, D_crit, Accentric factor, T_min, T_max, P_max, D_max
    #     self.eval(calculation, "1")

    #     # print(self.capture)
    #     res = re.findall(r"(\w*)\s* = \s*([-\d\.E]*)", self.capture)
    #     # print(res)
    #     if len(res) == 0:
    #         # print(self.capture) # show what trend output to console
    #         return None
    #     else:
    #         properties = {}
    #         for el in res:
    #             tmp_prop = Property(*el)
    #             properties[tmp_prop.name] = (tmp_prop.value, tmp_prop.unitname)
            
    #         if verbose:
    #             print(properties)
    #     return properties

    # def call_zz(self, calculation, verbose=False):
    #     raise(NotImplemented("Not currently required"))

    # def __init__(self, lib_path_str: Optional[str] = None):

    #     self._setup_library(lib_path_str)
    #     self._setup_library_functions()

    #     # init empty
    #     self.constants = None
    #     self.n_comp = 0
    #     self.fluid_name = None
    #     self.fluid_par = None
    #     self.is_fluid_name = None
    #     self.invalidate_density()
    #     self.invalidate_temperature()
    #     self.invalidate_pressure()
        
    #     self.get_eos_constants()

    # def __del__(self):
    #     """ this makes sure the dll allocated memory is freed after it is not needed"""
    #     self._deallocate_memory()
    #     # print('Destructor called') # DEBUG

    # def __str__(self):
    #     tmp_str =  f"EOS        : {self.eostype}\n"
    #     tmp_str += f"path       : {self.dll_path}\n"
    #     if self.is_fluid_name is True:
    #         tmp_str += f"fluid_name : {self.fluid_name.value}\n"
    #     elif self.is_fluid_name is False:
    #         tmp_str += f"fluid_par  : {self.fluid_par}\n"
    #     else:
    #         tmp_str += f"fluid      : {None}\n"
    #     tmp_str += f"n_comp     : {self.n_comp}\n"
    #     tmp_str += f"density    : {self.density}\n"
    #     tmp_str += f"temperature: {self.temperature}\n"
    #     tmp_str += f"pressure   : {self.pressure}\n"


    #     return tmp_str
    # # +++ SET METHODS +++
    # def set_n_comp(self, n_comp_new: int):
    #     """ set the new number of components if it differs from the used one"""

    # def set_fluid_by_name(self, fluid_name: List[str]):
    #     """ used to specify components by list of strings"""

    # def set_fluid(self, fluid_input: Optional[List[any]]):
    #     """ set one or multiple fluid names/parameters """

    # def set_fluid_by_par_single(self, fluid_parameters: Union[List[float], Tuple[float]], subst_index: int):
    #     """ set the single substances by their parameters"""

    # def set_temperature(self, temperature: Optional[float]):
    #     """ set the temperature [K] into the equation"""

    # def set_density(self, density: Optional[float]):
    #     """ set the molar density [mol/m^3] into the equation"""

    # def set_pfrac(self, pfrac, is_silent: bool = False):
    #     """ set the molar pfrac [#] into the equation"""

    # def set_dens_generalized(self, density: Optional[float], pfrac: Optional[float], is_silent: bool = False):
    #     """ unified way of setting density either regular or with packing fraction value"""

    # def set_pressure(self, pressure: Optional[float]):
    #     """ set the pressure [Pa] into the object equation does not operate with pressure input"""

    # def set_molar_ratio(self, molar_ratio: Optional[List[float]]):
     
    # # +++ INVALIDATE/UNSET METHODS +++
    # def invalidate_temperature(self):
    #     """ prevent accidental calculation with inpropper density setup by invalidating the density=None
    #     this makes special sence methods iterating over density """
    #     self.temperature = None

    # def invalidate_density(self):
    #     """ prevent accidental calculation with inpropper temperature setup by invalidating the temperature=None
    #     this makes special sence methods iterating over temperature """
    #     self.density = None

    # def invalidate_pressure(self):
    #     """ prevent accidental calculation with inpropper pressure setup by invalidating the pressure=None
    #     this makes special sence methods iterating over pressure """
    #     self.pressure = None

    # # +++ CALCULATE/GET DIRECT METHODS +++
    # def get_eos_constants(self) -> Dict[str, float]:
    #     """ get the dictionary with constant names and their corresponding values used in the eos"""
    #     if self.constants is None:
    #         r_gas = c_double()  # [J/(K*mol)]
    #         n_a = c_double()  # [1/mol]
    #         k_b = c_double()  # [J/K]
    #         pi = c_double()  # [1]

    #         self._get_constants(byref(r_gas), byref(n_a), byref(k_b), byref(pi))

    #         self.constants = {"gas_constant": r_gas.value,
    #                           "avogadro_constant": n_a.value,
    #                           "boltzmann_constant": k_b.value,
    #                           "pi": pi.value}
           
    #     return self.constants

    # def get_p(self, temperature: Optional[float] = None, density: Optional[float] = None, pfrac: Optional[float] = None,
    #           composition: Optional[List[float]] = None, substance: Optional[List[any]] = None,
    #           is_silent: bool = False):
    #     """ function that calculates pressure from temperature and density input"""

    # def get_mu_fug(self, temperature: Optional[float] = None, density: Optional[float] = None,
    #                composition: Optional[List[float]] = None, substance: Optional[List[any]] = None):


    # def get_h_s_g(self, temperature: Optional[float] = None, density: Optional[float] = None,
    #               composition: Optional[List[float]] = None, substance: Optional[List[any]] = None):


    # def get_d(self, temperature: Optional[float] = None, pressure: Optional[float] = None,
    #           substance: Optional[List[any]] = None, pfrac_estimates: Tuple[float, float, float] = (1e-10, 0, 0.4),
    #           is_silent: bool = False) -> Tuple[Optional[float], ...]:
    #     """ calculates density based on temperature and pressure, uses iterative Newton method over packing fraction
    #     does not need to converge all the time, try the density_estimates tuple to supply 0,1,2 root density estimate"""


    # def get_equil(self, temperature: Optional[float] = None, substance: Optional[List[any]] = None,
    #               pressure_estimate: float = 1e2, is_silent: bool = False, p_def_min: float = None) -> Tuple[
    #     Optional[float], Optional[float], Optional[float]]:
    #     """ return saturation properties or none in case of failure in form of (saturation_pressure, vapor_density, liquid_density)
    #      recomended pressure estimate is low number i.e for methane 100 Pa works fine as estimate"""
