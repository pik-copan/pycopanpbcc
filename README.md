# pycopanpbcc
Python scripts for analyzing Planetary Boundaries in a conceptual model of the Earth's Carbon Cycle including geoengineering by terrestrial carbon dioxide removal

# Associated publication
V. Heck, J.F. Donges, W. Lucht,
Collateral transgression of planetary boundaries due to climate engineering by terrestrial carbon dioxide removal, 
Earth System Dynamics 7, 783-796 (2016),
DOI: 10.5194/esd-7-783-2016.

Software DOI of pycopanpbcc release along with paper publication: 

[![DOI](https://www.zenodo.org/badge/72007144.svg)](https://www.zenodo.org/badge/latestdoi/72007144)

# General information
Python scripts run in Python 2.x and use standard Python packages numpy, scipy etc. 
Some scripts are designed to run on a parallel cluster system.

# Steps to reproduce results and figures 
1) The file scripts/main_model.py contains the main calibrated model (equations and parameters). If run as main file, Fig. 5 is produced and saved in home directory.

2) The file scripts/state_space_fig6.py calculates the different initial condition state space domains. If run as main file, Fig. 6 is produced and saved in home directory.

3) The file scripts/data_Fig7.py calculates the initial condition state space size for different domains for the variation of tCDR rate and threshold. If run as man file, data for the generation of Fig 7 are saved.

4) The file scripts/plots_fig7.py plots Fig7, based on data generated in script 3).

5) The files scripts/par_var RCP26.py, scripts/par_var RCP45.py, scripts/par_var RCP60.py, scripts/par_var RCP85.py generate data for Fig 8 for the respective RCP emission trajectories. They need to be run on a parallel cluster system and require the file mpi.py. 

6) The file scripts/plot_fig8.py generates Fig 8.

The data generated by scripts 3) and 5) are saved in /save/
