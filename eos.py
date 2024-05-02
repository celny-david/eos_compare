#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
perform a comparison of two eos under different settings
idea: is to verify the correctness of an eos against 
      another eos or another eos setting
{License_info}
"""

# Futures

# Built-in/Generic Imports
import os
import sys
from typing import List
# Libs
import click
from numpy import arange
# Own modules
from pc_saft.interface import pcsaft_eos
from TREND_4.interface import trend_eos
from TREND_4.conversion_table import is_valid_subst_name

__author__ = 'David Celny'
__copyright__ = 'Copyright {year}, {project_name}'
__license__ = '{license}'
__version__ = '{major}.{minor}.{rel}'
__maintainer__ = 'Dixiw'
__email__ = '{contact_email}'
__status__ = '{dev_status}' 



  
# pass_configuration = click.make_pass_decorator(class_instance_to_pass, ensure=True)
# NOTE example use
# @click.option('-c','--config_file', type=click.File('r'), default=None)
# @pass_configuration
# def sample(conf:initd.Configuration, config_file):

# ===== PARAMETERS ===== #
# name_to_eos = {
#                "trend": trend_eos,
#                "pc_saft": pcsaft_eos}
#NOTE dict of package: (eos,eos,...) entries
eos_available_packages={
    "trend": ("trend1","trend6"),
    "pcsaft": ("pcsaft",)
}

# ===== HELPERS ===== #
def parse_range(instr, output_type):
    """parse the range types of input """
    ndd = instr.count(":") #NOTE number of double dots
    if ndd !=2 or output_type is str:
        return [output_type(instr)]
    else:
        start,stop,step = instr.split(":")
        start = output_type(start)
        stop = output_type(stop)
        step = output_type(step)

        res_range = arange(start,stop,step)
        return res_range.tolist()

def parse_cs_args(ctx, param, value, output_type):    
    vals = str(value).split(',')
    # print(vals) # DEBUG
    if vals[0] == '':
        return []
    ret = []
    for val in vals:
        # print(type(val),val,output_type) #DEBUG
        ret += parse_range(val,output_type)
    return ret

def parse_cs_args_int(ctx, param, value):
    return parse_cs_args(ctx, param, value, int)

def parse_cs_args_float(ctx, param, value):
    return parse_cs_args(ctx, param, value, float)

def parse_cs_args_str(ctx, param, value):
    return parse_cs_args(ctx, param, value, str)

def get_eos(eos_name):
    # if eqname not in name_to_eos.keys():
    #     raise ValueError(f"EOS name: {eqname} not one from supported: {name_to_eos.keys()}")
    eos_name = eos_name.lower()
    if eos_name in ("trend", "trend1"):
        return trend_eos()
    if eos_name in ("trend6"):
        return trend_eos(eos_type=6)
    if eos_name in ("pcsaft", "pc_saft"):        
        pcsaft_path = os.path.join(os.getcwd(),"pc_saft","libpcsaft")
        return pcsaft_eos(lib_path_str=pcsaft_path)

    raise ValueError(f"EOS name: {eqname} not one from supported")

def print_data(data:List[float],format:str='16.6f'):
    """ utility for printing aligned data in compare"""
    for el in data:
        print(f"{el:{format}} ",end='')
    print() #NOTE to finish the line end

# ===== CLICK FUNCTIONALITY ===== #

@click.group()
def cli():    
    #NOTE place for general setup 
    pass

# === CLICK CLI GROUPS === #
@cli.command(short_help="help text here")
@click.argument("equation",          type=click.UNPROCESSED,
                callback=parse_cs_args_str)
@click.option("-s","--substance",    type=click.UNPROCESSED,             help="comma separated substance names",
              callback=parse_cs_args_str)
@click.option("-t","--temperature",  type=click.UNPROCESSED,             help="temperature input [K], csv or range start:stop:step",
              callback=parse_cs_args_float)
@click.option("-d","--density",      type=click.UNPROCESSED, default='', help="[opt] density input, [mol/m^3], csv or range start:stop:step",
              callback=parse_cs_args_float)
@click.option("-p","--pressure",     type=click.UNPROCESSED, default='', help="[opt] pressure input, [Pa], csv or range start:stop:step",
              callback=parse_cs_args_float)
@click.option("-o","--output",       type=click.UNPROCESSED,             help="comma separated output properties",
              callback=parse_cs_args_str)
@click.option("--diff","is_diff",    is_flag=True,                       help="show absolute difference vs first eos |1-2|")
@click.option("--rel" ,"is_rel",     is_flag=True,                       help="show relative difference vs first eos |2/1|")
@click.option("--format" ,"print_format",     type=str, default='16.6f', help="show relative difference vs first eos |2/1|")
def compare(equation, substance, temperature, density, pressure, output, is_diff, is_rel, print_format):
    """ perform comparison of two EoSs given by name (comma separated)"""
    # print(equation)
    # print(substance)
    # print(temperature)
    # print(density)
    # print(pressure)
    eoss = [get_eos(eqname) for eqname in equation]
    # for el in eoss: #DEBUG
    #     print(el)

    eos_padd = '10'

    if not is_valid_subst_name(substance, is_global=True):
        raise ValueError(f"Invalid substance {substance} not adhering to pc_saft naming convention.")


    print(f"Temp[K], Rho[mol/m^3], {'Eos_name':{eos_padd}}: {'Press[Pa]':>16} {'dP/dRho':>16} {'d2P/dRho2':>16}")    
    for t in temperature:
        for d in density:
            for ind,eos in enumerate(eoss):
                res = (eos.get_p(temperature=t,
                                 density=d,
                                 composition=[1.0,],
                                 # is_silent=True,
                                 substance=substance) )
                if ind == 0:
                    res1 = res
                    print(f"{t:7.1f}, {d:12.4f}, {eos.name_tag:{eos_padd}}: ",end='')
                else:
                    print(f"{' ':>22} {eos.name_tag:{eos_padd}}: ",end='')
                print_data(res, print_format)
                if is_diff and ind > 0:                    
                    diff = [abs(el1 - el2) for el1,el2 in zip(res1,res)]
                    print(f"{' ':>22} {'diff |1-'+str(ind+1)+'|':{eos_padd}}: ",end='')
                    print_data(diff, print_format)
                if is_rel and ind > 0:                    
                    rel = [el2/el1 for el1,el2 in zip(res1,res)]
                    print(f"{' ':>22} {'rel   1/'+str(ind+1):{eos_padd}}: ",end='')
                    print_data(rel, print_format)

@cli.command(short_help="return available EoS, EoS packages for comparison")
def list():
    align_factor = 12
    print(f"{'package name':<{align_factor}}: list of equations of state")
    print(f"{'-'*align_factor}:{'-'*27}")
    for package,eoss in eos_available_packages.items():
        print(f"{package:<{align_factor}}: ",end="")
        for eos in eoss:
            print(f"{eos}, ",end="")
        print()
    #DEBUG
    eoss = [get_eos(eqname) for eqname in equation]
    for el in eoss: 
        print(el)
    pass

@cli.group(help="for calculation from eos")
def calc():
    pass


# compare <what,what2> <compare params> <what to compare>

if __name__ == '__main__':
    cli()