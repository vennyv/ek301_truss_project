# EK301 Truss Design Project - MATLAB Code

**Boston University**  
**College of Engineering**  
**EK 301 Engineering Mechanics Spring 2024**  
**Truss Design Project**

## Project Overview

This repository contains the MATLAB code developed for the Truss Design Project in the course EK301 at Boston University, Spring 2024. The goal of this project is to design a truss capable of supporting a given load using specified materials. The design process is based on computational engineering analysis and optimization.

The key steps of the project include:
- Material characterization (buckling lab)
- Development of a truss analysis program
- Iterative design and optimization
- Final truss design and physical testing

### Project Specifications

- The truss must be a planar, simple truss and follow specific design constraints outlined in the project documentation.
- The goal is to minimize the total virtual cost while ensuring that the truss can support a minimum load of 32 oz.
- Truss member forces, reaction forces, and buckling behavior are computed using MATLAB.

### Key Dates
- **Buckling Lab**: March 1, 2024
- **Preliminary Design Report**: April 5, 2024
- **Final Design Report**: April 26, 2024
- **Truss Testing**: April 27, 2024

## Repository Structure

The MATLAB scripts included in this repository are used for truss analysis, uncertainty quantification, and optimization:

1. **EK301TrussProj.m**: The main script that computes the truss member forces, reaction forces, and total cost based on an input truss design. It solves the linear system for member forces using matrix methods.
  
2. **EK301TrussProjDesign1.m & EK301TrussProjDesign2.m**: Scripts that evaluate alternative truss designs and calculate the theoretical load-to-cost ratios. These scripts provide comparisons between different truss layouts and optimize for the best performance.
  
3. **EK301TrussProjPracticeProblem.m**: A practice script for testing and understanding truss analysis methods.
  
4. **uncertainty.m**: This script incorporates uncertainty in member buckling strengths using a fit formula based on experimental data. It calculates the failure load range and provides an estimated uncertainty in the truss's performance.

5. **EK301TrussProjectTest.m**: A test script to verify the accuracy of the truss analysis model by computing forces, reactions, and costs for a sample truss design.

## How to Run

1. Clone the repository:
   ```bash
   git clone https://github.com/your-repo/ek301_truss_project.git
   ```

2. Open MATLAB and navigate to the cloned directory.

3. Run the main analysis script for your truss design:
   ```matlab
   run('EK301TrussProj.m')
   ```

4. Modify the design matrices and vectors in the script to test different truss configurations. Example input matrices for various designs are provided in `EK301TrussProjDesign1.m` and `EK301TrussProjDesign2.m`.

## Example Output

The MATLAB scripts will output the following:
- Truss member forces (Tension/Compression)
- Reaction forces at supports
- Total cost of the truss design
- Predicted failure load (with uncertainty)

Example output from the test script:
```
Load: 32.0 oz
Member forces in oz:
m1: 45.32 (T)
m2: 12.15 (C)
...
Reaction forces in oz:
Sx1: 24.50
Sy1: 30.12
Sy2: 7.88
Cost of truss: $210.75
```

## Truss Design Process

The truss design is a structured process that includes:
1. **Material Characterization**: Data from the buckling lab is used to model the acrylic bar behavior under load.
2. **Model Development**: The MATLAB program predicts the performance of candidate truss designs.
3. **Design Iteration**: Multiple truss designs are evaluated to optimize the strength-to-cost ratio.
4. **Final Design & Testing**: The final design is constructed and tested for failure load against predictions.

## References

- [1] Gere, James M. *Mechanics of Materials*, 5th Edition.
- [2] Govindjee, Sanjay. *Engineering Mechanics of Deformable Solids*.
