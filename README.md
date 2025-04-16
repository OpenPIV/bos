# Background Oriented Schlieren for stratified liquid cases (Python Implementation)

## Overview

Background-oriented schlieren (BOS) is a technique for flow visualization of density gradients in fluids using the Gladstone–Dale relation between density and refractive index of the fluid.

BOS simplifies the visualization process by eliminating the need for expensive mirrors, lasers, and knife-edges. In its simplest form, BOS uses simple background patterns (randomly generated dot-pattern), an inexpensive strobe light source, and a high-speed digital camera.

This repository contains a Python implementation of the BOS technique for stratified liquid cases, translated from the original MATLAB code developed by Lilly Verso from Alex Liberzon's lab at Tel Aviv University.

## Installation

### Using Conda (Recommended)

```bash
# Create a new conda environment from the environment.yml file
conda env create -f environment.yml

# Activate the environment
conda activate bos_env
```

### Using pip

```bash
# Install required packages from requirements.txt
pip install -r requirements.txt
```

## Running the Demo

To run the BOS demo, use the following command:

```bash
python openpiv_bos_demo.py
```

This will process the example data and generate visualizations of the results.

## Code Structure

### Main Files
- `openpiv_bos_demo.py` - Demo script for BOS processing
- `BOS_MAIN_ROUTINE.py` - Main script for BOS processing
- `bos_poisson_solver.py` - Solves the Poisson equation for BOS
- `bos_remapping.py` - Remaps images according to displacement field
- `bos_correlation_openpiv.py` - Performs PIV correlation between images
- `imwarp.py` - Warps image with flow field
- `parameters.py` - Parameters for BOS processing
- `gladstone_dale.py` - Applies Gladstone-Dale relation
- `poisson_direct_mod.py` - Alternative Poisson solver using direct method
- `create_grid.py`, `create_rhs.py`, `crop_field.py`, etc. - Utility functions

### OpenPIV Python Package
- `openpiv_python/` - Contains the Python implementation of the OpenPIV algorithms

### Poisson Test Python Package
- `poisson_test_python/` - Contains the Python implementation of the Poisson solvers

## Testing

Test scripts are provided to verify the Python implementation:

```bash
# Test the imwarp function
python test_imwarp.py

# Test various BOS functions
python test_bos_translation.py

# Test the OpenPIV Python functions
python test_openpiv_python.py

# Test the Poisson test Python functions
python test_poisson_test_python.py
```

## Requirements

- Python 3.6+
- NumPy
- SciPy
- Matplotlib
- Pillow (PIL)
- OpenCV (cv2)

## How to cite this work

Verso, L. and Liberzon, A. "Background Oriented Schlieren in a Density Stratified Fluid", Rev. Sci. Instrum. 86, 103705 (2015)
 
```bibtex
@article{Verso:2015,  
   author = "Verso, Lilly and Liberzon, Alex",  
   title = "Background oriented schlieren in a density stratified fluid",  
   journal = "Review of Scientific Instruments",  
   year = "2015",  
   volume = "86",  
   number = "10",   
   eid = 103705,  
   url = "http://scitation.aip.org/content/aip/journal/rsi/86/10/10.1063/1.4934576",    
   doi = "http://dx.doi.org/10.1063/1.4934576" 
}
```
