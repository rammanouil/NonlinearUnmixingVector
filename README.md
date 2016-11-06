# GLUP

This demo reproduces the results reported in the paper

>"GLUP: Yet another algorithm for blind unmixing of hyperspectral data," <br />
> R. Ammanouil, A. Ferrari, C. Richard, and D. Mary, <br />
>6th Workshop on Hyperspectral Image and Signal Processing: Evolution in Remote Sensing (WHISPERS), Lausanne, pp. 1-4 (2014)


1. Synthetic Data GLUP example : This folder allows to reproduce the first simulation showing the image of the estimated abundance matrix. It contains:
  * Function CreateData.m that allows to create synthetic data
  * Matlab variable endmembers.mat contains 8 endmembers from USGS library
  * Matlab variable 200pxl8end40dB1.mat an example of a synthetic data created using CreateData.m
  * Function GLUP.m that performs GLUP
  * Script main.n that loads variable, set parameters for GLUP and calls GLUP

2. Synthetic Data Compare Approaches : This folder allows to reproduce the second simulation in the paper,
    where we compare the performance of GLUP, NFINDR and SDSOMP by repeating the simulation 100 times. It contains :
  * Matlab variable endmembers.mat
  * Function GLUP.m
  * Function SDSOMP.m that searches for endmembers based on a Greedy approach in reference "Greedy algorithms
      for pure pixels identification in hyperspectral"
  * Function ZMD_NFINDR.m based on reference ": An algorithm for fast autonomous spectral endmember determination
      in hyperspectral data"
  * Script main2.m that loads endmembers, creates synthetic data, repeats 100 times the three algorithms and notes
      the percentage of accurately detected endmembers. The algorithm notes the results in "GNS_test1.txt"
  * Text document "GNS_test1.txt" contains the simulation results reported in table 1 in paper

3. Real Data Compare Approaches : This folder allows to reproduce the third simulation done using a real HSI, namely,
    Pavia Center image. It contains :
  * Function Angle_Difference.m computes spectral angle between a matrix of column wise data and its estimated version
  * Function Compute_RMSE.m computes RMSE between a matrix of column wise data and its estimated version
  * Function GLUP.m
  * Script main.m which loads Pavia Center image, then performs GLUP, NFINDR and SDSOMP, and shows the abundance
      matrix of each methods endmembers. the algorithm also computes RMSE and avg spectral angle.
  * Function qpas.m to perform FCLS
  * Function SDSOMP.m
  * Function ZMD_NFINDR
  * Matlab variable Pavia.mat

