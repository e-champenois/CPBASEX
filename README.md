CPBASEX
=======

pBASEX Abel Transform Inversion without Polar Rebinning... in MATLAB and Python

Some features:
-implementation of the pBASEX Abel inversion algorithm by sampling on a cartesian grid (i.e. straight from an image) rather than rebinning into polar coordinates and introducing some error.
-use of singular value decomposition to speed up the least-squares fitting step.
-basis function regularized (also fast) or pixel weighted (slightly slower) least-squares also possible.
-parallelized and efficient Abel transformed basis function calculations.

Inversion code for MATLAB and Python is available in the CpBASEX folders, while the example directory has examples of how to calculate the integrals needed for the inversion (saveGData) or invert some sample images and plot the results (basicAnalysis).

A full write-up on Abel inversion, pBASEX in general, details of this particular implementation, and comparison with other methods will be available soon...