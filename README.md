[![Build Status](https://travis-ci.org/frank-y-liu/SPRNT.svg?branch=master)](https://travis-ci.org/frank-y-liu/SPRNT)
[![License](https://img.shields.io/badge/License-EPL%201.0-red.svg)](https://opensource.org/licenses/EPL-1.0)
<img alt="Docker Pulls" src="https://img.shields.io/docker/pulls/frankliu1/sprnt">
### SPRNT 

=====

SPRNT is a river dynamics simulation software package. Technically it is a fully dynamic
Saint Venant Equation solver. Quite a few nice techniques are developed to make SPRNT not 
only fast but also capable of handling large networks. The code is mainly written by
Dr. Frank Liu at IBM Research, as the result of academic collaborations with Professor Ben
Hodges in the Civil Engineering Department of University of Texas at Austin. It is released
under Eclipse Public License. IBM owns the copyright.

## Content: 
* src/        directory contains the core SPRNT routine source code
* spt/        contains the code of a front-end to run SPRNT as command-line 
* examples/   contains several netlist examples
* demo/       contains an example on how to run SPRNT through API
* doc/        contains a user's manual on the syntax of the netlist
* ThirdParty/ contains modules to build solver library as dynamic-linked library  

## Install
See INSTALL on how to build and install the software

## Docker Image
Now you can run SPRNT in container. See README.docker for details. Follow this link to pull the image from dockerhub: https://hub.docker.com/r/frankliu1/sprnt.

## Contact
For bugs/questions/comments/critiques, please send email through GitHub.

## References

* Liu, Frank, and Ben R. Hodges. "[Applying microprocessor analysis methods to river
network modelling](https://dx.doi.org/10.1016/j.envsoft.2013.09.013)." Environmental Modelling & Software 52 (2014): 234-252.

* Hodges, Ben R., and Frank Liu. "[Rivers and Electric Networks: Crossing Disciplines in
Modeling and Simulation](https://dx.doi.org/10.1561/1000000033)" Foundations and Trends in Electronic Design Automation 8.1
(2014): 1-116.

* Yu, Cheng-Wei, Frank Liu and Ben R. Hodges. "[Consistent initial conditions for the Saint-Venant equations in river network modeling](https://doi.org/10.5194/hess-21-4959-2017)" Hydrology and Earth System Science 21.9(2017) 4959.  

* Yu, Cheng-Wei, Ben R. Hodges and Frank Liu, "[A new form of the Saint-Venant equations for variable topography](https://hess.copernicus.org/articles/24/4001/2020/)" Hydrology and Earth System Science 24.8(2020) 4001.
