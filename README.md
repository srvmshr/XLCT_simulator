### Introduction

Simulation for selective excitation based X-ray Luminescent CT.  


### Acknowledgements 

G Pratx, Stanford University

#### Essential Requirements
* MATLAB 2008a+
* MATLAB Parallel Computing Toolbox

Note: This repo was tested on MATLAB R2011a. There might have been breaking changes 
in MATLAB's API. Hence, any such loss of funtionality is to be taken into consideration 
before running these codes. The main code running the simulation is `test_run.m`. This is 
a non-GPU optimized code, but parallelization was possible on a cluster with 64 CPU nodes.
Try running this on a cluster with `qsub`.

#### File Details

| Contents 		| Description 							|
|----------- 		|-------------							|
|**compinter.m** 	| Perfoms computation of grid intersection points.		|
|**genlines.m**  	| Generates lines of predefined spacing for XLCT excitation	|
|**gen_phan.m**		| Generates the phantom						| 
|**mc321.c**		| Tiny Monte-carlo simulator					|
|**raytrace.m**		| Raytrace the projection					|
|**raytrace_recon**	| Backprojection raytracing					|
|**recon.m**		| Reconstruction function					|
|**runsim.m**		| Runs sinogram simulation					|
|**test_run.m**		| Main program							|
|**sens.png**		| Calculated Sensitivity map					|
|**water** 		| Water background (Monte Carlo)				|
