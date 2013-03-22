comparingGSalgorithms
=====================
MATLAB codes for comparing geostatistical simulation algorithms.
 Please do the following to use the codes: 
1. Copy-paste the codes in a specific directory of your computer. 
2. Go to ...\ ComparingGSAlgorithms\ComparingGSAlgorithms folder and open test.m file. 
3. uncomment and run the code for each case. 

9 important functions are listed below, explaining how this algorithm works in detail.
Functions: 
  1.  DistMtrx = calculateModelVar_MPH(realizations,TI,pyramid)
  
  This function is used to calculate distance matrix between realizations, as well as between realizations and training image in 2D binary case. MPH approach is used in this function.
  First input:  “realizations” is all the realizations needed: e.g.101*101*50
  Second input: “TI” is the training image
  Third input:  “pyramid” is the level of Pyramid of realizations and training image
  Output:       ”DistMtrx”: distance matrix e.g. 10*51*51 (10 is the level of Pyramid)

   Notice: sometimes it might be out of contiguous memory, so rebooting or restarting MATLAB may work. When possible, running the code in 64-bit windows will solve such a problem.



  2. DistMtrx = calculate3DModelVar_MPH(realizations,TI,pyramid)
  This function is used to calculate distance matrix between realizations, as well as between realizations and training image in 3D binary case. MPH approach is used in this function.

  First input: “realizations” is all the realizations needed: e.g.69*69*39*50
  Second input: “TI” is the training image: e.g. 69*69*39
  Third input:  “pyramid” is the level of Pyramid of realizations and training image
  Output:       “DistMtrx”: distance matrix e.g. 10*51*51 (10 is the level of Pyramid)

   Notice: sometimes it might be out of contiguous memory, so rebooting or restarting MATLAB may work. When possible, running the code in 64-bit windows will solve such a problem.


  3. DistMtrx = calculateModelVar_CHP(realizations,TI,tempSize);
   This function is used to calculate distance matrix between realizations, as well as between realizations and training image in continuous case. CHP approach is used in this function.

   First input: “realizations” is all the realizations needed: e.g.101*101*50
   Second input: “TI” is the training image
   Third input: “tempSize” is the optimal size of template for training image
   Output:      “DistMtrx”: e.g. distance matrix e.g. 10*51*51 (10 is the level of Pyramid)


  4. show2DModelVar2(DisMtrx1,DisMtrx2,TI,realizations1,realizations2)
   This function is used to show the variability of the two sets of 2D realizations at different resolutions in the form of MDS plots.

   First input: DisMtrx1 is the distance matrix for the first set of realizations
   Second input: DisMtrx2 is the distance matrix for the second sets of realizations
  Third input: “TI” is the training image e.g. 101*101
  Forth input: “realizations1” is the first set of realizations: e.g. 101*101*50
  Fifth input: “realizations2” is the second set of realizations

  5. show2DModelVar3(DisMtrx1,DisMtrx2,DisMtrx3,TI,realizations1,realizations2,realizations3)
   This function is used to show the variability of the three sets of 2D realizations at different resolutions in the form of MDS plots.

   First input: DisMtrx1 is the distance matrix for the first set of realizations
   Second input: DisMtrx2 is the distance matrix for the second sets of realizations
   Third input: DisMtrx3 is the distance matrix for the third sets of realizations
  Forth input: “TI” is the training image e.g. 101*101
  Fifth input: “realizations1” is the first set of realizations: e.g. 101*101*50
  Sixth input: “realizations2” is the second set of realizations
  Seventh input: “realizations3” is the third set of realizations

  6. show3DModelVar2(DisMtrx1,DisMtrx2,TI,realizations1,realizations2)
   This function is used to show the variability of the two sets of 3D realizations at different resolutions in the form of MDS plots.

   First input: DisMtrx1 is the distance matrix for the first set of realizations
   Second input: DisMtrx2 is the distance matrix for the second sets of realizations
  Third input: “TI” is the training image e.g. 69*69*39
  Forth input: “realizations1” is the first set of realizations: e.g. 69*69*39*50
  Fifth input: “realizations2” is the second set of realizations

  7. show3DModelVar3(DisMtrx1,DisMtrx2,DisMtrx3,TI,realizations1,realizations2,realizations3)
   This function is used to show the variability of the three sets of 3D realizations at different resolutions in the form of MDS plots.

   First input: DisMtrx1 is the distance matrix for the first set of realizations
   Second input: DisMtrx2 is the distance matrix for the second sets of realizations
   Third input: DisMtrx3 is the distance matrix for the third sets of realizations
  Forth input: “TI” is the training image e.g. 69*69*39
  Fifth input: “realizations1” is the first set of realizations: e.g. 69*69*39*50
  Sixth input: “realizations2” is the second set of realizations
  Seventh input: “realizations3” is the third set of realizations


  8. RatioMtrx=QuantifyDist2(DisMtrx1,DisMtrx2)
   This function is used to calculate the ratio of between-realization distance and within-realization distance according to two different distance matrix. 
   First input: DisMtrx1 is the distance matrix for the first set of realizations
   Second input: DisMtrx2 is the distance matrix for the second sets of realizations
   Output: RatioMtrx is a 3x2 matrix to represent ratios.
                            

  9. RatioMtrx=QuantifyDist3(DisMtrx1,DisMtrx2,DisMtrx3)
   This function is used to calculate the ratio of between-realization distance and within-realization distance according to three different distance matrix. 
   First input: DisMtrx1 is the distance matrix for the first set of realizations
   Second input: DisMtrx2 is the distance matrix for the second sets of realizations
   Third input: DisMtrx3 is the distance matrix for the second sets of realizations
   Output: RatioMtrx is a 3x3 matrix to represent ratios
