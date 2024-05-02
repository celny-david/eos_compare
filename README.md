# eos_compare
Perform comparison of calculation from two Equation of state (EoS) for purpose of verification of predicitive capabilities of the EoS's. Beware that this varsion is intended as showcase. For operation the required libraries for TREND 4.1 and PC-SAFT are required.

## requirements
As this is research related tool originaly intended for personal use during development of PC-SAFT EoS very specific requirements have to be fulfilled for comparison to work.
### TREND
Because of cross-compatibility issues with older version of TREND the interface to TREND 4.1 is constructed with issuing calls into main program using console cpature to gather the results. In this way foreigh codebase could be incorporated without modification of the external codebase.
#### to make it work
The path to the main program then needs to be specified in `TREND_4/interface.py` in variable `programpath` on line `44`.
### PC-SAFT
Because this tool was developed alongside PC-SAFT a conveniet interface is included making comunication much faster and more reliable. For the required libpcsaft library and parameter_list.txt one can refer to [this PC-SAFT](https://github.com/celny-david/pc_saft_eos) project.
#### to make it work
There are placeholder files that needs to be replaced with files obtained from compliation of the PC-SAFT dynamic_lib target and `parameter_list.txt` found in `res` folder


## Running eos_compare
Here is result of help showing available options.
``` shell
Usage: eos.py compare [OPTIONS] EQUATION

  perform comparison of two EoSs given by name (comma separated)

Options:
  -s, --substance TEXT    comma separated substance names
  -t, --temperature TEXT  temperature input [K], csv or range start:stop:step
  -d, --density TEXT      [opt] density input, [mol/m^3], csv or range
                          start:stop:step
  -p, --pressure TEXT     [opt] pressure input, [Pa], csv or range
                          start:stop:step
  -o, --output TEXT       comma separated output properties
  --diff                  show absolute difference vs first eos |1-2|
  --rel                   show relative difference vs first eos |2/1|
  --format TEXT           show relative difference vs first eos |2/1|
  --help                  Show this message and exit.
```

### supported equations:
 - trend/trend1 = TREND multiparametric helmholtz type of EoS
 - trend6 = TREND PCP-SAFT type of EoS
 - pcsaft = own fortran implementation of PCP-SAFT EoS


### example call
 how to call the comparison for two eqations from pc-saft and pc-saft implemented in TREND for `methane` within temperature range `<200,300>` spaced with `20K` for single molar density `0.1mol/m^3`
``` shell
	./eos.py compare pcsaft,trend6 -s methane -t 200.0:300.0:20.0 -d 0.1
```