# RRMSD metric

## Goal
We want to find optimal charges to parametrize the SMAMP molecules with. For this, we use HORTON charge fitting, which has two free parameters: rhoref and sigma.
The HORTON algorithm minimizes the weighted error between the eletrostatic potential (ESP) calculated via DFT calculation, and the ESP obtained from point charges.
Rhoref determines the electron density at which the most weight will be assigned.
Sigma is the broadness of the weight, i.e. how fast the weights decay away from rhoref.
In the original paper, values of $10**{-9}$ (rhoref) and 0.8 (sigma) were used, and it was suggested that this would generalize to other molecules

## Problem
We have seen that the default values do not behave well for the charged structures Sarah wants to simulate. This project aims to find reasonable values.

## Solution
To evaluate which parameters are better suited, we do some parameter sweeps.
For the evaluation of the sweeps, we need a function mapping the charges to real numbers.
In the original paper, the authors used RRMSD (Relative Root Mean Square Deviation) between DFT-ESP and Pointcharge-ESP.
Sarah wants to use RRMS between constrained and unconstrained charges instead. This metric is defined as:
$RRMSD = \sqrt{\frac{\sum_i^N{(q^{con}_i - q^{uncon}_i)^2}}{N}}}$

, with a sum over all atoms/residue i, with a positive value in units of elementary charge [e].
This metric puts more weight on outliers than a non-quadratic sum would.

Evaluating the RRMSD over the 2D-space of (rhoref,sigma), we can evaluate how close our chosen parameter values are to the minimum.
