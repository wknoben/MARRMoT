# MARRMoT
Modular Assessment of Rainfall-Runoff Models Toolbox - Matlab code for 46 conceptual hydrologic models.

<p align="center">
<img src="Figures/logo.jpg" alt="MARRMoT logo" width="200"/>
</p>

MARRMoT is a novel rainfall-runoff model comparison framework that allows objective comparison between different conceptual hydrological model structures. 
The framework provides Matlab code for 46 unique model structures, standardized parameter ranges across all model structures and robust numerical implementation of each model.
The framework is provided with extensive documentation, a User Manual and several workflow scripts that give examples of how to use the framework.
MARRMoT is based around individual flux functions and aggregated model functions, allowing a wide range of possible applications.

If you have any questions about using or running the code, or are willing to contribute, please contact wouter.knoben[-at-]usask.ca. 

[![DOI](https://zenodo.org/badge/161804123.svg)](https://zenodo.org/badge/latestdoi/161804123)

## Changes after peer review
MARRMoT v1.2 has been accepted through peer review. Since then, users have found various bugs which are corrected on the current master branch. Summary:

- 'interflow_9' was missing a non-negativity constraint. This has been added and included in MARRMoT v1.3.
- Water balance calculations did not properly account for time step sizes different than 1 day. This has been corrected on the master branch but not been released yet.
- Models m05, m15, m37 and m44 did not properly account for time step size in certain flux calculations. This has been corrected on the master branch but not been released yet.
- Workflow_example_4 now works under Octave, after Octave update 5.2.0 and code contribution by Mustafa Kemal Türkeri. This has been integrated on the master branch but not been released yet.
- A new model m47 (m_47_IHM19_16p_4s) has been added after a contribution by Clara Brandes and her supervisors. 
- Several additional efficiency metrics have been added, thanks to Thomas Whöling.
- Efficiency metrics now accept a warmup period (number of initial time steps to ignore when calculating efficiency metrics) as an optional argument. Should be backwards compatible with existing scripts. Thanks to Thomas Whöling.

## Getting Started
These instructions will help you install a copy of MARRMoT and run a few example cases. 

### Requirements
MARRMoT was developed on Matlab version 9.2.0.538062 (R2017a) and the Optimization Toolbox is required (tested with version 7.6). 
If using Octave, MARRMoT was tested on Octave 4.4.1 and requires the 'optim' package.  
The following instructions assume that MARRMoT will be used with Matlab. 

Note that the function `circshift()` that is used by routing routines has markedly different behaviour in Matlab 2016b and higher compared to previous versions. Routing results will be unreliable in Matlab 2016a and below but will **_not_** generate any warnings or error messages. User discretion is advised.

### Install
Download a copy of the files from this repository (note: do not use the folder 'Octave' when using Matlab) and extract the files in an appropriate directory.

### Try an example application
- Open Matlab
- Add the folder 'MARRMoT' and its subfolders 'Functions', 'Models' and 'User Manual' to the Matlab path
- **Note:** Ensure that the folder "Octave" is **not** included, and remove this folder if it is
- Navigate Matlab's current folder to './MARRMoT/User Manual'
- Open the script 'workflow_example_1.m'
- Run the script by pressing F5
- Repeat with 'workflow_example_2.m' and 'workflow_example_3.m' ('workflow_example_4.m' shows a calibration example and takes a bit longer)

The User Manual provides further details.

## Documentation
The article describing MARRMoT development will soon be submitted to the scientific journal 'Geoscientific Model Development'.
This paper, its Supporting Material and the User Manual cover the following topics:

- **Paper**: rationale behind MARRMoT development, best practices used during development, summary of included model structures and an example application of all structures to simulate streamflow in a single catchment. https://doi.org/10.5194/gmd-12-2463-2019
- **Supporting Material (full)**: https://www.geosci-model-dev.net/12/2463/2019/gmd-12-2463-2019-supplement.pdf
- **Supporting Material (section 2)**: detailed description of each model structure, giving Ordinary Differential Equations and constitutive functions for each model store
- **Supporting Material (section 3)**: translation of constitutive functions (fluxes) to Matlab code
- **Supporting Material (section 4)**: overview of Unit Hydrograph code
- **Supporting Material (section 5)**: rationale behind generalised parameter ranges (use of these ranges is optional)
- **User Manual**: covers a variety of topics including (i) understanding model files, (ii) application examples, (iii) creating a new model or flux function, and (iv) Octave-specific instructions. https://github.com/wknoben/MARRMoT/blob/master/MARRMoT%20User%20manual%20v1.2.pdf

## Model structure summary
MARRMoT model structures are based on a wide variety of different models. 
However, do to the standardised format of this framework, MARRMoT models resemble, but are not the same as the models they are based on.
In addition to a range of unnamed models, the following models provided inspiration for MARRMoT:

- FLEX-Topo
- IHACRES
- GR4J
- TOPMODEL
- SIMHYD
- VIC
- LASCAM
- TCM
- TANK
- XINANJIANG
- HYMOD
- SACRAMENTO
- MODHYDROLOG
- HBV-96
- MCRM
- SMAR
- NAM
- HYCYMODEL
- GSM-SOCONT
- ECHO
- PRMS
- CLASSIC

## License
MARRMoT is licensed under the GNU GPL v3 license - see the LICENSE file for details.

## DOIs of previous releases
- v1.2: dx.doi.org/10.5281/zenodo.3235664
- v1.1: dx.doi.org/10.5281/zenodo.2677728
- v1.0: dx.doi.org/10.5281/zenodo.2482542 

## Acknowledgements
MARRMoT could not have been made without the effort that many hydrologists have put into development of their models. Their effors are gratefully acknowledged. Special thanks are extended to:
- Philip Kraft for finding a bug in the flux smoothing code during peer review; 
- Sebastian Gnann for suggesting various quality of life fixes; 
- Clara Brandes for finding and suggesting a fix for a bug in the water balance calculations and implementing m47;
- Koen Jansen for suggesting various improvements and correcting parameter descriptions;
- Mustafa Kemal Türkeri for making workflow_example_4 operational in Octave; and for performing extensive testing of MARRMoT in Matlab and Octave;
- Thomas Whöling for suggesting various additional efficiency metrics and a possible implementation for warmup periods.
