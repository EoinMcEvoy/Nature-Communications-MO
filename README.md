# Nature-Communications-MO
Code from McEvoy et al. (2020) Nature Communications

*****************************************************
This repository contains three files:

1. A MATLAB (last tested on ver 2019a) file to simulate mechano-osmotic behavior of two connected cells
"two_cell_dynamic_R.m"

2. A MATLAB (last tested on ver 2019a) file to simulate mechano-electro-osmotic behavior of two connected cells
"two_cell_dynamic_MEO.m" 

3. A COMSOL Multiphysics (last tested on ver5.4a) model to simulate a cell cluster behavior
"Cell_cluster_R.mph"

* All files and software were tested in Windows 10.
** Code commented within "two_cell_dynamic_R.m" for detailed explanation of analysis.

*****************************************************
Matlab code demo:
1. Install MATLAB (https://www.mathworks.com/products/matlab.html). The installation usually takes approximately 1 hour.
2. Open two_cell_dynamic.m and click "Run".
* New windows will then open with the simulated results

two_cell_dynamic_R.m:
This file is configured to reproduce simulated mechano-osmotic behavior of connected cells with operative gap junctions [Figures 2B-D and S1]
For GJ inhibition [Figures 2E-G and S2] : In the default file set wg = 0 [line 5] and Lpg = 0 [line 6].
* Figure colors are set to correspond with Figure 2B-D. 
* axes may need to be adjusted to correspond with manuscript figures

two_cell_dynamic_MEO.m:
This file is configured to reproduce simulated mechano-electro-osmotic behavior of connected cells with operative gap junctions 
Steady state predictions will reflect those reported in Figure S6. 
For GJ inhibition (or single cells) [Figure S5] : In the default file set w = 0 [line 23] and Lpg = 0 [line 24].
* axes may need to be adjusted to correspond with manuscript figures

*****************************************************
COMSOL simulation demo:
1. Install COMSOL Multiphysics (https://www.comsol.com/). The installation usually takes approximately 1.5 hours.
2. Navigate to "Model Builder" [left panel], click on "Study" and then click on "Parametric sweep". 
3. Navigate to "Settings" [middle panel] and click on "Compute".
* All computational results reported within Figure 3 will be visible within the "Results" tab in the "Model Builder" panel.

Note: This file is configured to reproduce results shown in Figure 3.
All other results may be reproduced my modifying parameters as in accordance with the manuscript. Briefly:
- Figure 4 (GJ inhibitor) : Navigate to the "Parametric Sweep" settings and set wg = 0 and Lpg = 0. Click "Compute".
- Figure 4 (stress release) : Navigate to the "Parametric Sweep" settings and set sig_g_max = 400. Click "Compute".

*****************************************************
* To ensure accurate reproduction of results we recommend you reopen the unedited files provided before modifying parameters and settings.
* The demos may take 1-5 minutes to run depending on the machine. 

*****************************************************
For any questions regarding these files, please contact Eoin McEvoy or Vivek B. Shenoy.
