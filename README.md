# Background Oriented Schlieren for stratified liquid cases


### Background Oriented Schlieren

From Wikipedia, the free encyclopedia
Background-oriented schlieren (BOS) is a novel technique for flow visualization of density gradients in fluids using the Gladstone–Dale relation between density and refractive index of the fluid.

BOS simplifies the visualization process by eliminating the need for the use of expensive mirrors, lasers and knife-edges. In its simplest form, BOS makes use of simple background patterns of the form of a randomly generated dot-pattern, an inexpensive strobe light source and a high speed digital camera.

In its initial stages of implementation, it is mostly being used as a qualitative visualization method. Further research into developing the technique will enable quantitative visualization of fluid flows, in various applications including cryogenic flows, supersonic and hypersonic flows, biomedical device flow visualization etc.

Read more on <https://en.wikipedia.org/wiki/Background-oriented_schlieren_technique>


### What is so special about liquid cases 



### What is so special about stratified cases


### Who wrote the software

The software is written mainly by Lilly Verso from Alex Liberzon's lab at Tel Aviv University, using 
OpenPIV package as a platform. The software utilizes the OpenPIV Matlab package for the cross-correlation analysis 
(essentially a stripped version of PIV analysis) and OpenPIV - pressure package for Poisson solver ideas. Additional 
Poisson solvers were tested, using public domain Matlab codes. The authors of these packages are gratefully acknowledged.

### Dependencies

1. OpenPIV-Matlab http://www.openpiv.net/openpiv-matlab/
2. OpenPIV-Pressure http://www.openpiv.net/openpiv-pressure/
3. IMWRAP - <https://github.com/animesh-garg/videoSeg/blob/master/src/classic_NL_Code/utils/imwarp.m>  under <https://github.com/animesh-garg/videoSeg> 
4. 2D Poisson equation - <http://www.mathworks.com/matlabcentral/fileexchange/38090-2d-poisson-equation/content/Poisson_equation_2D.m>


### How to cite this work: 

Verso, L. and Liberzon, A. "Background Oriented Schlieren in a Density Stratified Fluid", <http://arxiv.org/abs/1506.08889>