# MARRMoT
Modular Assessment of Rainfall-Runoff Models Toolbox - Matlab code for 46 conceptual hydrologic models.

<p align="center">
<img src="Figures/logo.jpg" alt="MARRMoT logo" width="200"/>
</p>

MARRMoT is a novel rainfall-runoff model comparison framework that allows objective comparison between different conceptual hydrological model structures. 
The framework provides Matlab code for 46 unique model structures, standardized parameter ranges across all model structures and robust numerical implementation of each model.
The framework is provided with extensive documentation, a User Manual and several workflow scripts that give examples of how to use the framework.
MARRMoT is based around individual flux functions and aggregated model functions, allowing a wide range of possible applications.

If you have any questions about using or running the code, or are willing to contribute, please contact w.j.m.knoben[-at-]bristol.ac.uk.

## Getting Started
These instructions will help you install a copy of MARRMoT and run a few example cases. 

### Requirements
MARRMoT was developed on Matlab version 9.2.0.538062 (R2017a) and the Optimization Toolbox is required (tested with version 7.6). 
If using Octave, MARRMoT was tested on Octave 4.4.1 and requires the 'optim' package. 
The following instructions assume that MARRMoT will be used with Matlab. 

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

- **Paper**: rationale behind MARRMoT development, best practices used during development, summary of included model structures and an example application of all structures to simulate streamflow in a single catchment
- **Supporting Material 2**: detailed description of each model structure, giving Ordinary Differential Equations and constitutive functions for each model store
- **Supporting Material 3**: translation of constitutive functions (fluxes) to Matlab code
- **Supporting Material 4**: overview of Unit Hydrograph code
- **Supporting Material 5**: rationale behind generalised parameter ranges (use of these ranges is optional)
- **User Manual**: covers a variety of topics including (i) understanding model files, (ii) application examples, (iii) creating a new model or flux function, and (iv) Octave-specific instructions

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

## Acknowledgements
MARRMoT could not have been made without the effort that many hydrologists have put into development of their models. Their effors are gratefully acknowledged.
