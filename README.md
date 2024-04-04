# mine_thermal_energy_transport
Evaluate the techno-economic feasibility of shallow geothermal energy transport; see XXXXX (pending link for article).

## Requirements
### Overall 
 - Python 3.10 (code has not been tested with other Python versions).
 - Python package dependencies:
 -   Numpy
 -   Pandas
 -   Scipy
 -   Matplotlib
 -   Optional: [TSAM](https://github.com/FZJ-IEK3-VSA/tsam) for time series aggregation.
 
### Model Execution
 - Python package dependencies:
   - [COMANDO](https://jugit.fz-juelich.de/iek-10/public/optimization/comando) and its dependencies (Optimisation framework).
   - [Psweep](https://github.com/elcorto/psweep/tree/0.9.0) 0.10.0 or lower (Parametric analysis and housekeeping).
 - [GUROBI](https://www.gurobi.com/) 10.0 or higher

## Setup
 - Check for Python version and package dependencies (virtual environment or separated Conda environment is recommended)
 - SETUP GUROBI solver installation, further details can be found in the COMANDO [readthedocs](https://comando.readthedocs.io/en/latest/interfaces.html#interfaces)
  
---
## Required Inputs

- Thermal load time series. The original study considered time series for single-family houses generated using [RC_BuildingSimulator](https://github.com/architecture-building-systems/RC_BuildingSimulator)
- Ambient temperature time series. The original study used TMY data from [PV-GIS](https://re.jrc.ec.europa.eu/pvg_tools/en/)
- Ground temperature time series (for pipe losses calculation). The original study derived the ground temperature from the ambient temperature, as indicated in the manuscript.
- Horizon of analysis and interest rate (n = 30 years and i = 0.03 by default).
- Cost of components and commodities.


## Usage
### Base_Code
The code includes a sample case for a decentralised thermal system, considering conditions for the Rheinish Coal Mining Area (Rheinisches Braunkohlerevier), NRW, Germany region and a non-insulated PE transmission pipeline. The code is executed by running the study.py file. All the required inputs and functions are in the _data and _script folders. Results will be stored in a newly generated _results folder, and a summary of the results will be stored in a pandas data frame inside a pickle file (.pk).

Further component examples are found in the [COMANDO](https://jugit.fz-juelich.de/iek-10/public/optimization/comando) list of examples.

## Considerations
- COMANDO is solver agnostic; minor modifications to the code should be needed to run the code with other solvers, although the authors haven't tested this so far. Part of the code used for post-processing relies on GUROBI outputs, which would need to be adjusted as well.
The individual heat pumps/RCAC components are reversible. If time aggregation is used, be careful with simultaneous heating and cooling demand, as this will result in an unfeasible model.
- The study file is set up to run using Python multiprocessing capabilities. If this isn't desired, the code must be changed accordingly.
