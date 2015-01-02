CPBASEX
=======

[MATLAB] PBASEX Abel Inversion without Polar Rebinning

This README needs some work...

Instructions:

(1) Save all of the files in a directory belonging to your MATLAB path.

(2) Try to run saveGData. It might error because the two functions it calls findG and findGinv use parfor loops instead of for loops to make speed up the computations by doing them in parallel. If they do error because of parfor not being available, open findG.m and findGInv.m and replace the word "parfor" in the code with "for", and run saveGData again. This will take a long time (if you want to test it faster, change the first 3 lines to:

gData.x = 0:63;
gData.y = x;
gData.k = 0.5:2:62.5;

If this worked, there should be a new .mat file in your current directory with all the needed data for inversions.

(3) 
