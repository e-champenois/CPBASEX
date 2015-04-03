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

(3) Load up your raw image(s) to invert and add up either all or some of the quadrants to be left with a 2-D (or 3-D if multiple images) array where the first entry is the center of the image (upper-left of the remaining quadrant). Also, the size of the image should correspond to the size of gData.x and gData.y above, so either chop what you need from the image, resize the image, or redo the gData calculation with a gData.x and gData.y that suit your image.

(4) Run the inversion: call the pbasex function as done in its example in its comment section.
