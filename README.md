# VC-BayesianEstimation

[![License](https://img.shields.io/badge/license-BSD%203--clause-green.svg)](LICENSE)
[![GitHub release](https://img.shields.io/badge/version-v1.5.0-blue.svg)](https://github.com/vcurdia/VC-BayesianEstimation/releases/tag/v1.5.0)

Codes used to estimate a Dynamic Stochastic General Equilibrium (DSGE) model
using Bayesian Estimation techniques.


# Requirements

## Matlab (R)
The codes were tested using Matlab (R) R2012a with the following toolboxes
- Symbolic Toolbox
- Statistical Toolbox
- Optimization Toolbox

## LaTeX
LaTeX is used by some tools to compile certain documents.

`epstopdf`, included in most LaTeX releases, is used by some tools.

## Additional codes needed
- [VC-Tools](https://github.com/vcurdia/VC-Tools)
  by
  [Vasco Cúrdia](http://www.frbsf.org/economic-research/economists/vasco-curdia/), 
  version 
  [v1.5.0](https://github.com/vcurdia/VC-Tools/releases/tag/v1.5.0)
- [gensys](http://sims.princeton.edu/yftp/gensys/)
  by [Chris Sims](http://www.princeton.edu/~sims/)
- [optimize](http://dge.repec.org/codes/sims/optimize/)
  by [Chris Sims](http://www.princeton.edu/~sims/)
- [KF](http://sims.princeton.edu/yftp/Times09/KFmatlab/)
  by [Chris Sims](http://www.princeton.edu/~sims/)


## Usage example

The script `SetDSGE.m` is an example of how to setup the model and estimate it
using this package.

The main structure for setup
1. Set file names for the data input and the output
2. Set parameters list and priors
3. Set list of observation variables
4. Set list of State space variables
5. Set list of iid shocks
6. Generate symbolic variables
7. Construct any necessary auxiliary definitions (optional)
8. Set observation equations
9. Set state Equations

The above 9 steps will set up the model. After setup the Bayesian estimation
proceeds by finding the mode of the posterior using `MaxPost.m` and then
generating MCMC samples, using `MCMC.m`. Analysis of estimation results is done
with `MCMCAnalysis.m`.

See the example `SetDSGE.m` for basic options and how to call the sequence of
steps in more detail. 

Reports in pdf format are generated along the way.


# Additional Information

Each of the functions and scripts contains help at the beginning of the codes,
including options and flags. This help can be accessed in an interactive Matlab
(R) session using the `help` or `doc` commands: 
```
help SetDSGE
``` 
or 
```
doc SetDSGE
```

