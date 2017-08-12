# Notes on different files

## Debye Calculation
[debye.py](debye.py):
parallelizes Debye calculation over DCD frames

[multiJob_debye.py](multiJob_debye.py):
same as [debye.py](debye.py) but used with multi PDB/DCD pairs

[parallel_q_debye.py](parallel_q_debye.py):
Basically the same as debye.py except parallaizes over Q rather than frames. Useful for large PDB files with only one frame.

## Full I(Q) Calculation

[scatteringCalc.py](scatteringCalc.py):
Full exponential calculator. Same inputs as [debye.py](debye.py).

## Cube Subtraction

[scatteringCalc-Results.ipynb](scatteringCalc-Results.ipynb):
- `doFit`: does the cube fitting
- `cubeScat`: calls the cube scattering method from sas-models (ripped out of SASView)

## Radial Distribution Function, g(r)

[GOFR.ipynb](GOFR.ipynb):
calculate g(r) from S(q)

## Create PDBs

[createPDBs.ipynb](createPDBs.ipynb):
Create a single PDB by repositioning several copies of a single PDB.


# Notes from [Ian Hunt-Isaak](https://github.com/ianhi)
## Requires

nbstripout
`pip install nbstripout`

Do not use git kracken when commiting .ipynb files as it doesn't seem to run the proper nbstripout command.
Rather you should commit from the command line.

## Data Files
Are located on Ian's computer and should be placed in `/home/data/sascaclc_pbc`
