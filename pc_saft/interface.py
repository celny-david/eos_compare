# -*- coding: utf-8 -*-

"""
This is the interface for the PC-SAFT dynamic library call.
Expose the so/dll functionality for two supported systems {LINUX,WIN}.
"""
__author__ = 'David Celný'
__version__ = '0.0.1'

# Built-in/Generic Imports
import os
import sys
import platform
# […]
from ctypes import *
from typing import Optional, Tuple, List, Dict, Callable, Union

# Libs
import numpy as np
from scipy.optimize import newton, root_scalar, minimize, bisect

# Own modules
from eos_abc import Eos
try:
    import pc_saft.root_iter as ri
except:    
    import root_iter as ri

def is_iterable(object):
    """ verify if the object is of iterable type"""
    try:
        iter(object)
        return True 
    except TypeError:
        return False

class pcsaft_eos(Eos):
    lib_path_default = "libpcsaft"  # default name for the library
    num_of_substace_par = 11
    eostype = "PC-SAFT"
    name_tag = "pcsaft"

    def _add_lib_ext(self, lib_path: str, platform: str) -> str:
        _, ext = os.path.splitext(lib_path)
        if ext == '':
            if platform == "Linux":
                lib_path += '.so'
            elif platform == "Windows":
                lib_path += '.dll'
            else:
                raise NotImplemented(f"Platform {platform} is not supported!")
        return lib_path

    def _setup_library(self, lib_path_str: Optional[str]):
        """ group the setting/opening of library on both platforms"""
        self.is_deallocate_memory = False
        self.ps = platform.system()
        if lib_path_str is None:
            dll_name = self._add_lib_ext(self.lib_path_default, self.ps)
        else:
            dll_name = self._add_lib_ext(lib_path_str, self.ps)
        if self.ps == "Linux":
            dll_path = os.path.realpath(dll_name) # or .dll or .so
            if os.path.exists(dll_path):
                self.lib = cdll.LoadLibrary(dll_path)
                self.dll_path = dll_path
            else:
                self.is_deallocate_memory = False
                raise FileNotFoundError(f"library {dll_path} not found, check if the file specified exists.")

        elif self.ps == "Windows":
            dll_path = os.path.realpath(dll_name) # or .dll or .so
            if os.path.exists(dll_path):
                self.lib = CDLL(dll_path, winmode=0)
                self.dll_path = dll_path
                # libtest = WinDLL(dll_path, winmode=0)
            else:
                self.is_deallocate_memory = False
                raise FileNotFoundError(f"library {dll_path} not found, check if the file specified exists.")
        # print("dll_path", dll_path) #DEBUG
        self._setup_library_functions()

    def _setup_library_functions(self):
        """ set up the library functions and their arguments"""

        # LIBRARY FUNCTIONS and types
        self._get_constants = self.lib.get_constants
        self._get_constants.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double)]

        self._deallocate_memory = self.lib.deallocate_memory
        self._deallocate_memory.argtypes = []
        
        self._parse_input = self.lib.parse_input
        self._parse_input.argtypes = [c_char_p, c_int]

        self._set_n_comp = self.lib.set_n_comp
        self._set_n_comp.argtypes = [c_int]

        self._set_component_name = self.lib.set_component_name
        self._set_component_name.argtypes = [c_char_p, c_int]

        self._set_component_par = self.lib.set_component_param
        self._set_component_par.argtypes = [c_int, POINTER(c_double)]

        self._set_density = self.lib.set_density
        self._set_density.argtypes = [c_double]
        
        self._set_pfrac = self.lib.set_pfrac
        self._set_pfrac.argtypes = [c_double, POINTER(c_double)]

        self._set_temperature = self.lib.set_temperature
        self._set_temperature.argtypes = [c_double]

        self._set_molar_ratio = self.lib.set_molar_ratio
        self._set_molar_ratio.argtypes = [POINTER(c_double)]

        self._press_calc = self.lib.press_calc
        self._press_calc.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double)]

        self._chem_pot_calc = self.lib.chem_pot_calc
        self._chem_pot_calc.argtypes = [POINTER(c_double), POINTER(c_double)]

        self._enthalpy_entropy_calc = self.lib.enthalpy_entropy_calc
        self._enthalpy_entropy_calc.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double)]

        self._all_calc = self.lib.all_calc
        self._all_calc.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double),
                                   POINTER(c_double), POINTER(c_double),
                                   POINTER(c_double), POINTER(c_double), POINTER(c_double)]
        self.is_deallocate_memory = True

    def __init__(self, lib_path_str: Optional[str] = None):

        self._setup_library(lib_path_str)

        # init empty
        self.constants = None
        self.n_comp = 0
        self.fluid_name = None
        self.fluid_par = None
        self.is_fluid_name = None
        self.invalidate_density()
        self.invalidate_temperature()
        self.invalidate_pressure()
        
        self.get_eos_constants()

    def __del__(self):
        """ this makes sure the dll allocated memory is freed after it is not needed"""
        #NOTE to make sure this is tried only when it makes sense (i.e. not in case when library is not found)
        if self.is_deallocate_memory:
            self._deallocate_memory()
        # print('Destructor called') # DEBUG

    def __str__(self):
        tmp_str =  f"EOS        : {self.eostype}\n"
        tmp_str += f"path       : {self.dll_path}\n"
        methods = [method for method in dir(pcsaft_eos) if callable(getattr(pcsaft_eos, method)) and not method.startswith("_")]
        tmp_str += f"methods    : {methods},\n"

        if self.is_fluid_name is True:
            tmp_str += f"fluid_name : {self.fluid_name.value}\n"
        elif self.is_fluid_name is False:
            tmp_str += f"fluid_par  : {self.fluid_par}\n"
        else:
            tmp_str += f"fluid      : {None}\n"
        tmp_str += f"n_comp     : {self.n_comp}\n"
        tmp_str += f"density    : {self.density}\n"
        tmp_str += f"temperature: {self.temperature}\n"
        tmp_str += f"pressure   : {self.pressure}\n"

        return tmp_str

    def __eq__(self, other):
        if (self.lib_path_str != other.lib_path_str):
            return False
        if (self.eostype != other.eostype):
            return False
        if (self.fluid_name != other.fluid_name):
            return False
        return True
    # +++ SET METHODS +++
    def set_n_comp(self, n_comp_new: int):
        """ set the new number of components if it differs from the used one"""
        if self.n_comp != n_comp_new:
            self._set_n_comp(n_comp_new)
            self.n_comp = n_comp_new

    def set_fluid_by_name(self, fluid_name: List[str]):
        """ used to specify components by list of strings"""
        
        fluid_str = ''
        for name in fluid_name:
            fluid_str += name
            fluid_str += ','
        fluid_str = fluid_str[:-1]

        self.fluid_name = c_char_p(fluid_str.encode('utf-8'))

        n_comp_new = len(fluid_name)
        self.set_n_comp(n_comp_new)
        self._set_component_name(self.fluid_name, len(self.fluid_name.value))
        self.is_fluid_name = True
        # print(self.fluid_name.value.split(b','))

    def set_fluid(self, fluid_input: Optional[List[any]]):
        """ set one or multiple fluid names/parameters """
        if fluid_input is None:
            if self.is_fluid_name is None:
                raise ValueError("Warning: No fluid specified or set -> set_fluid function before call or as argument")
            if self.is_fluid_name is True and self.fluid_name is None:
                raise ValueError(
                    "Warning: No fluid name specified or set -> set_fluid([fluid_name]) function before call or as argument")
            if self.is_fluid_name is False and self.fluid_par is None:
                raise ValueError(
                    "Warning: No fluid parameters specified or set -> set_fluid([fluid_parameters]) function before call or as argument")
        else:
            if len(fluid_input) <= 0:
                raise ValueError("The passed list has to be nonempty!")

            if type(fluid_input[0]) is str:
                self.set_fluid_by_name(fluid_input)

            elif type(fluid_input[0]) is str:
                self.set_n_comp(1)
                self.set_fluid_by_par_single(fluid_input, 0)
                self.fluid_par = fluid_input
                self.is_fluid_name = False

            elif type(fluid_input[0]) is list or type(fluid_input[0]) is tuple:
                n_comp_new = len(fluid_input)
                self.set_n_comp(n_comp_new)
                self.fluid_par = []
                for par_ind, par in enumerate(fluid_input):
                    self.set_fluid_by_par_single(par, par_ind)
                    self.fluid_par.append(par)
                self.is_fluid_name = False
            else:
                raise ValueError(f"Unexpected zero element type: {type(fluid_input[0])} in the list")
    def set_fluid_by_par_single(self, fluid_parameters: Union[List[float], Tuple[float]], subst_index: int):
        """ set the single substances by their parameters"""
        if len(fluid_parameters) != self.num_of_substace_par:
            raise ValueError(
                f"Wrong number of substance {subst_index} parameters {len(fluid_parameters)} != {self.num_of_substace_par}")
        fluid_parameters_np = np.asarray(fluid_parameters, dtype=c_double)
        # print(subst_index, fluid_parameters)
        self._set_component_par(subst_index + 1, fluid_parameters_np.ctypes.data_as(POINTER(c_double)))

    def set_temperature(self, temperature: Optional[float]):
        """ set the temperature [K] into the equation"""
        if temperature is None:
            if self.temperature is None:
                raise ValueError(
                    "Warning: No value specified or set -> temperature= argument or set_temperature function before call")
        else:
            if is_iterable(temperature):
                temperature = float(temperature[0])
            self.temperature = temperature

            self._set_temperature(c_double(self.temperature))

    def set_density(self, density: Optional[float]):
        """ set the molar density [mol/m^3] into the equation"""
        if density is None:
            if self.density is None:
                raise ValueError(
                    "Warning: No value specified or set -> density= argument or set_density function before call")
        else:
            if is_iterable(density):
                density = float(density[0])
            if density <= 0.0:
                self.density = 1e-2
                print(f"Warning: density <=0 : setting new value: {self.density}")
            else:
                self.density = density

            self._set_density(c_double(self.density))

    def set_pfrac(self, pfrac, is_silent: bool = False):
        """ set the molar pfrac [#] into the equation"""
        if pfrac is None:
            raise ValueError("Error: No value pfrac is")
        else:
            # raise NotImplemented("not yet implemented")
            dens_out = c_double()
            if is_iterable(pfrac):
                pfrac = float(pfrac[0])

            if pfrac <= 0.0:
                pfrac = 1e-16
                if is_silent is False:
                    print(f"Warning: unreasonably low pfrac <=0 is changed to: {pfrac}")

            closest_packing = np.pi / (3 * np.sqrt(2))
            if pfrac > closest_packing:
                pfrac = closest_packing
                if is_silent is False:
                    print(f"Warning: unreasonably high pfrac > closesct pcaking is changed to: {pfrac}")
                
            self._set_pfrac(c_double(pfrac), byref(dens_out))
            if dens_out.value < 0.0:
                raise RuntimeError(f"Failed to set packing fraction with error code: {dens_out.value}")

            return dens_out.value

    def set_dens_generalized(self, density: Optional[float], pfrac: Optional[float], is_silent: bool = False):
        """ unified way of setting density either regular or with packing fraction value"""
        if pfrac is not None:
            if density is not None:
                raise ValueError("Error unable to accept simultaneous input of density and pfrac.")
            density = self.set_pfrac(pfrac, is_silent)
            if density > 0.0:
                self.density = density
            else:
                raise ValueError("Invalid negative value of density")
        else:
            self.set_density(density)

    def set_pressure(self, pressure: Optional[float]):
        """ set the pressure [Pa] into the object equation does not operate with pressure input"""
        if pressure is None:
            if self.pressure is None:
                raise ValueError(
                    "Warning: No value specified or set -> pressure= argument or set_pressure function before call")
        else:
            if is_iterable(pressure):
                pressure = float(pressure[0])
            self.pressure = pressure


    def set_molar_ratio(self, molar_ratio: Optional[List[float]]):

        if molar_ratio is None:
            if self.molar_ratio is None:
                raise ValueError(
                    "Warning: No value(s) specified or set -> molar_ratio= argument or set_molar_ratio function before call")
        else:
            self.molar_ratio = np.asarray(molar_ratio, dtype=c_double)

            self._set_molar_ratio(self.molar_ratio.ctypes.data_as(POINTER(c_double)))
        
    # +++ INVALIDATE/UNSET METHODS +++
    def invalidate_temperature(self):
        """ prevent accidental calculation with inpropper density setup by invalidating the density=None
        this makes special sence methods iterating over density """
        self.temperature = None

    def invalidate_density(self):
        """ prevent accidental calculation with inpropper temperature setup by invalidating the temperature=None
        this makes special sence methods iterating over temperature """
        self.density = None

    def invalidate_pressure(self):
        """ prevent accidental calculation with inpropper pressure setup by invalidating the pressure=None
        this makes special sence methods iterating over pressure """
        self.pressure = None

    # +++ CALCULATE/GET DIRECT METHODS +++
    def get_eos_constants(self) -> Dict[str, float]:
        """ get the dictionary with constant names and their corresponding values used in the eos"""
        if self.constants is None:
            r_gas = c_double()  # [J/(K*mol)]
            n_a = c_double()  # [1/mol]
            k_b = c_double()  # [J/K]
            pi = c_double()  # [1]

            self._get_constants(byref(r_gas), byref(n_a), byref(k_b), byref(pi))

            self.constants = {"gas_constant": r_gas.value,
                              "avogadro_constant": n_a.value,
                              "boltzmann_constant": k_b.value,
                              "pi": pi.value}
           
        return self.constants

    def get_p(self, temperature: Optional[float] = None, density: Optional[float] = None, pfrac: Optional[float] = None,
              composition: Optional[List[float]] = None, substance: Optional[List[any]] = None,
              is_silent: bool = False):
        """ function that calculates pressure from temperature and density input, and density derivatives
            p [Pa], pp [Pa*m**3/mol], ppp [Pa*m**6/mol**2]"""
        p = c_double()
        pp = c_double()
        ppp = c_double()

        self.set_fluid(substance)
        self.set_molar_ratio(composition)
        self.set_temperature(temperature)
        self.set_dens_generalized(density=density, pfrac=pfrac, is_silent=is_silent)

        #NOTE disable console output from the following function
        if is_silent:
            old_stdout = sys.stdout # backup current stdout
            print("Nothing in between this start")
            sys.stdout = open(os.devnull, "w")

        self._press_calc(byref(p), byref(pp), byref(ppp))

        if is_silent:
            sys.stdout = old_stdout
            print("Nothing in between this stop")

        return (p.value, pp.value, ppp.value)

    def get_mu_fug(self, temperature: Optional[float] = None, density: Optional[float] = None,
                   composition: Optional[List[float]] = None, substance: Optional[List[any]] = None):

        self.set_fluid(substance)
        self.set_molar_ratio(composition)
        self.set_temperature(temperature)
        self.set_density(density)
        mu = np.empty((self.n_comp),dtype=c_double())
        fug = np.empty((self.n_comp),dtype=c_double())

        self._chem_pot_calc(mu.ctypes.data_as(POINTER(c_double)),
                            fug.ctypes.data_as(POINTER(c_double)))

        # print(fug) # DEBUG
        if self.n_comp == 1:
            mu = mu[0]
            fug = fug[0]
        return (mu, fug)

    def get_h_s_g(self, temperature: Optional[float] = None, density: Optional[float] = None,
                  composition: Optional[List[float]] = None, substance: Optional[List[any]] = None):

        h = c_double()
        s = c_double()
        g = c_double()

        self.set_fluid(substance)
        self.set_molar_ratio(composition)
        self.set_temperature(temperature)
        self.set_density(density)

        self._enthalpy_entropy_calc(byref(h), byref(s), byref(g))

        return (h.value, s.value, g.value)

    # +++ SECOND LEVEL CALCUALTE METHODS +++
    # using the direct interface methods
    # def get_d(self, temperature:Optional[float]=None, pressure:Optional[float]=None, substance:Optional[List[any]] = None, density_estimates:Tuple[float,float,float] = (10,10,10), is_silent:bool = False) -> Tuple[Optional[float], ...]:
    #     """ calculates density based on temperature and pressure, uses iterative Newton method over density
    #     does not need to converge all the time, try the density_estimates tuple to supply 0,1,2 root density estimate"""

    #     self.set_fluid(substance)
    #     self.set_temperature(temperature)
    #     self.set_pressure(pressure)
        
    #     if self.n_comp != 1:
    #         raise NotImplemented(f"get_d is implemented for pure fluid only and 1!=n_comp: {self.n_comp}")
                
    #     # NOTE - Ideal gas estimate of density if required
    #     # dens_est = self.pressure/(self.constants["gas_constant"]*self.temperature)
                
    #     # print(f"density estimate: {density_estimates}") #DEBUG
    #     # TODO remove pressure from self
    #     roots = ri.d_root(self.get_p, target_pressure = self.pressure)
        
    #     roots.is_silent = is_silent

    #     # === ROOT 0 -optimally is vapor ===
    #     d_root0 = roots.get_root0(density_estimates[0])
    #     # === ROOT 1 -meta/un stable in all cases===
    #     d_root1 = roots.get_root1(density_estimates[1], d_root0)
    #     if d_root1 is None:
    #         return (d_root0, None, None)
    #     # === ROOT 2 -optimally is liquid ===
    #     self.invalidate_density()
    #     d_root2 = roots.get_root2(density_estimates[2], d_root0, d_root1)
    #     if d_root2 is None:
    #         if d_root1 > d_root0:
    #             return (d_root0, d_root1, None)
    #         elif d_root1 < d_root0:
    #             return (d_root1, d_root0, None)
    #         else:
    #             return (d_root1, None, None)
    #     else:
    #         return sorted((d_root0, d_root1, d_root2))

    def get_d(self, temperature: Optional[float] = None, pressure: Optional[float] = None,
              substance: Optional[List[any]] = None, pfrac_estimates: Tuple[float, float, float] = (1e-10, 0, 0.4),
              is_silent: bool = False) -> Tuple[Optional[float], ...]:
        """ calculates density based on temperature and pressure, uses iterative Newton method over packing fraction
        does not need to converge all the time, try the density_estimates tuple to supply 0,1,2 root density estimate"""

        self.set_fluid(substance)
        self.set_temperature(temperature)
        self.set_pressure(pressure)
        
        if self.n_comp != 1:
            raise NotImplemented(f"get_d is implemented for pure fluid only and 1!=n_comp: {self.n_comp}")
                
        # NOTE - Ideal gas estimate of density if required
        # dens_est = self.pressure/(self.constants["gas_constant"]*self.temperature)
                
        # print(f"density estimate: {density_estimates}") #DEBUG
        # TODO remove pressure from self
        # roots = ri.d_root(self.get_p, target_pressure = self.pressure)
        # roots.is_silent = is_silent

        roots = ri.pfrac_root(self.get_p, pfrac_to_dens=self.set_pfrac, target_pressure=self.pressure)

        roots.is_silent = is_silent

        # === ROOT 0 -optimally is vapor ===
        d_root0 = roots.get_root0(pfrac_estimates[0])
        # === ROOT 1 -meta/un stable in all cases===
        d_root1 = roots.get_root1(pfrac_estimates[1])
        # === ROOT 2 -optimally is liquid ===
        d_root2 = roots.get_root2(pfrac_estimates[2])
        
        # print(d_root0, d_root1, d_root2) # DEBUG

        if (d_root0 is not None) and (d_root2 is not None) and (abs(d_root0 - d_root2) < 1e-12):
            # print(f"set root to None: {abs(d_root0 - d_root2)}") # DEBUG
            d_root0 = None

        return (d_root0, d_root1, d_root2)


    def get_equil(self, temperature: Optional[float] = None, substance: Optional[List[any]] = None,
                  pressure_estimate: float = 1e2, is_silent: bool = False, p_def_min: float = None) -> Tuple[
        Optional[float], Optional[float], Optional[float]]:
        """ return saturation properties or none in case of failure in form of (saturation_pressure, vapor_density, liquid_density)
         recomended pressure estimate is low number i.e for methane 100 Pa works fine as estimate"""

        self.set_fluid(substance)
        self.set_temperature(temperature)
        if self.n_comp != 1:
            raise NotImplemented(f"get_equil is implemented for pure fluid only and 1!=n_comp: {self.n_comp}")
      
        self.last_it = [0,0,0]

        # NOTE function that backtrack itself when it overshoots 
        def f(p_est):
            d1 = None
            d3 = None
            cnt = 0
            maxcnt = 10
            while (d1 is None or d3 is None) and cnt < maxcnt:
                p_used = float(p_est - cnt * (p_est - self.last_it[0]) / maxcnt) # NOTE consecutive maxcnt equidistant return
                # p_used = self.last_it[0] + (p_est - self.last_it[0])/ 2**cnt # NOTE consecutive halving
                
                if p_def_min is not None and p_used <= 0 :
                    if is_silent is False:
                        print(f"Forced correction of negative pressure: {p_used} to: {p_def_min*10**cnt}")
                    p_used = p_def_min*10**cnt
                d1, d2, d3 = self.get_d(temperature=self.temperature, pressure=p_used, is_silent=True)
                # print(d1,d3)
                if cnt > 0 and is_silent is False:
                    print(f"For T={self.temperature} return from: {p_est} by {cnt}/{maxcnt} to: {p_used}")
                cnt += 1

            # if p_def_min is not None and p_used <= 0 :
            #     if is_silent is False:
            #         print(f"Forced correction of negative pressure: {p_used} to: {p_def_min}")
            #     p_used = p_def_min
            #     d1, d2, d3 = self.get_d(temperature=self.temperature, pressure=p_used, is_silent=True)    

            # print(f"p: {p_used}")
            # print(f"{d1}, {d2}, {d3}")
            # print(f"{d1:.20}, {d2:.20}, {d3:.20}")

            if d1 is None or d3 is None:                
                raise RuntimeError(
                    f"For T={self.temperature}K, P_est={pressure_estimate}, iterative algorithm ended in unphysical state")
            # if d1 is None:
            #     mu_vap = 0.0
            #     fug_vap = 0.0
            # else:
            mu_vap, fug_vap = self.get_mu_fug(temperature=self.temperature, density=d1)
            mu_liq, fug_liq = self.get_mu_fug(temperature=self.temperature, density=d3)
            self.last_it[0] = p_used
            self.last_it[1] = d1
            self.last_it[2] = d3
            # DEBUG SECTION
            # print("vap:",fug_vap, mu_vap)
            # print("liq:",fug_liq, mu_liq)
            # print("val:",fug_liq - fug_vap, mu_liq - mu_vap)
            # print("div:",fug_vap/fug_liq -1.0, mu_vap/mu_liq - 1.0)
            # print("div-:",1.0 - fug_vap/fug_liq, 1.0 - mu_vap/mu_liq)
            # print()
            # various value forms
            # return (fug_liq - fug_vap)
            # return (1.0 - fug_vap/fug_liq)
            return (fug_vap/fug_liq - 1.0)
            # return (mu_liq - mu_vap)
            # return (mu_liq/mu_vap - 1.0)

        try:
            # result = minimize(f,
            #                   x0=pressure_estimate,
            #                   options={'maxiter':100},
            #                   method='COBYLA',
            #                   # method='Nelder-Mead',
            #                   tol=10 ** (-1)
            #                   )
            # if result.success is True:
            #     print(result.x)
            #     return round(result.x[0], ndigits=1), self.last_it[1], self.last_it[2]

            solution, report = newton(f,
                                      x0=pressure_estimate,
                                      # maxiter=100,
                                      full_output=True,
                                      tol=10 ** (-1)
                                      )
            #print(report) # DEBUG
            if report.converged is True:
                return round(solution, ndigits=1), self.last_it[1], self.last_it[2]
            else:
                return (None, None, None)
        except RuntimeError as e:
            print(e)
            return (None, None, None)
            
    def get_equil_rf(self, temperature: Optional[float] = None, substance: Optional[List[any]] = None,
                  pressure_estimates: float = [1.0,5e7], is_silent: bool = False, p_def_min: float = None) -> Tuple[
        Optional[float], Optional[float], Optional[float]]:
        """ return saturation properties or none in case of failure in form of (saturation_pressure, vapor_density, liquid_density)
         recomended pressure estimate is low number i.e for methane 100 Pa works fine as estimate"""

        self.set_fluid(substance)
        self.set_temperature(temperature)

        if self.n_comp != 1:
            raise NotImplemented(f"get_equil is implemented for pure fluid only and 1!=n_comp: {self.n_comp}")
      
        self.last_it = [0,0,0]

        # NOTE function that backtrack itself when it overshoots 
        def f(p_est):

            p_used = float(p_est)

            d1, _, d3 = self.get_d(temperature=self.temperature, pressure=p_used, is_silent=True)

            if d1 is None and d3 is None:
                raise RuntimeError(f"For T={self.temperature}K, P_est={p_used}, iterative algorithm ended in unphysical state")
            
            if d1 is None:
                mu_vap = 0.0 
                fug_vap = 0.0
                mu_liq, fug_liq = self.get_mu_fug(temperature=self.temperature, density=d3)
                mu_liq *= -1
            elif d3 is None:
                mu_vap, fug_vap = self.get_mu_fug(temperature=self.temperature, density=d1)
                mu_liq, fug_liq = (0.0, 1.0)
            else:
                mu_vap, fug_vap = self.get_mu_fug(temperature=self.temperature, density=d1)
                mu_liq, fug_liq = self.get_mu_fug(temperature=self.temperature, density=d3)

            # d1 = None
            # d3 = None
            # cnt = 0
            # maxcnt = 10
            # while (d1 is None or d3 is None) and cnt < maxcnt:
            #     p_used = float(p_est - cnt * (p_est - self.last_it[0]) / maxcnt) # NOTE consecutive maxcnt equidistant return
            #     # p_used = self.last_it[0] + (p_est - self.last_it[0])/ 2**cnt # NOTE consecutive halving
                
            #     if p_def_min is not None and p_used <= 0 :
            #         if is_silent is False:
            #             print(f"Forced correction of negative pressure: {p_used} to: {p_def_min*10**cnt}")
            #         p_used = p_def_min*10**cnt
            #     d1, d2, d3 = self.get_d(temperature=self.temperature, pressure=p_used, is_silent=True)
            #     # print(d1,d3)
            #     if cnt > 0 and is_silent is False:
            #         print(f"For T={self.temperature} return from: {p_est} by {cnt}/{maxcnt} to: {p_used}")
            #     cnt += 1
            
            # # print(f"p: {p_used}")
            # # print(f"{d1}, {d2}, {d3}")
            # # print(f"{d1:.20}, {d2:.20}, {d3:.20}")

            # if d1 is None or d3 is None:                
            #     raise RuntimeError(
            #         f"For T={self.temperature}K, P_est={pressure_estimate}, iterative algorithm ended in unphysical state")
            # # if d1 is None:
            # #     mu_vap = 0.0
            # #     fug_vap = 0.0
            # # else:
            # mu_vap, fug_vap = self.get_mu_fug(temperature=self.temperature, density=d1)
            # mu_liq, fug_liq = self.get_mu_fug(temperature=self.temperature, density=d3)
            self.last_it[0] = p_used
            self.last_it[1] = d1
            self.last_it[2] = d3
            # DEBUG SECTION
            # print("vap:",fug_vap, mu_vap)
            # print("liq:",fug_liq, mu_liq)
            # print("val:",fug_liq - fug_vap, mu_liq - mu_vap)
            # print("div:",fug_vap/fug_liq -1.0, mu_vap/mu_liq - 1.0)
            # print("div-:",1.0 - fug_vap/fug_liq, 1.0 - mu_vap/mu_liq)
            # print()
            # various value forms
            # return (fug_liq - fug_vap)
            # return (1.0 - fug_vap/fug_liq)
            # return (fug_vap/fug_liq - 1.0)
            # return (mu_liq - mu_vap)
            return (mu_vap - mu_liq)
            # return (mu_liq/mu_vap - 1.0)

        try:
            x, result = bisect(f,
                            a = pressure_estimates[0],
                            b = pressure_estimates[1],
                            full_output = True)

            if result.converged is True:
                # print(result.root)
                return round(result.root, ndigits=1), self.last_it[1], self.last_it[2]

            # solution, report = newton(f,
            #                           x0=pressure_estimate,
            #                           # maxiter=100,
            #                           full_output=True,
            #                           tol=10 ** (-1)
            #                           )
            # #print(report) # DEBUG
            # if report.converged is True:
            #     return round(solution, ndigits=1), self.last_it[1], self.last_it[2]
            # else:
            #     return (None, None, None)
        except RuntimeError as e:
            print(e)
            return (None, None, None)