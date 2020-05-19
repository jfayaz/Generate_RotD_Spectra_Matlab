# Generate RotD Spectra using Matlab


This code develops the RotD50 Sa and RotD100 Sa Spectra of the Bi-Directional  Ground Motion records provided in the 'GMs.mat' file which must be in the current folder. 

GMs.mat must contain 3 variables

    'acc1'              --> n x 1 Cell Structure containing GM history in Direction 1 

    'acc2'              --> n x 1 Cell Structure containing GM history in Direction 2

    'dt'                --> n x 1 array containing dts for the GM histories


INPUT:

This codes provides the option to have 3 different regions of developing the Spectra of ground motions with different period intervals (discretizations)
 
The following inputs within the code are required:
 
     'Int_T_Reg_1'        --> Period Interval for the first region of the Spectrum 

     'End_T_Reg_1'        --> Last Period of the first region of the Spectrum (where to end the first region)

     'Int_T_Reg_2'        --> Period Interval for the second region of the Spectrum 

     'End_T_Reg_2'        --> Last Period of the second region of the Spectrum (where to end the second region)

     'Int_T_Reg_3'        --> Period Interval for the third region of the Spectrum 

     'End_T_Reg_3'        --> Last Period of the third region of the Spectrum (where to end the third region)

     'Damp'               --> Value of Damping for SDOF system

     'Plot_Spectra'       --> whether to plot the generated Spectra of the ground motions (options: 'Yes', 'No')    

 
OUTPUT:

The output will be provided in cell structure 'SPECTRA.mat' file which contain RotD50 and RotD100 spectra of the GMs in the same order as given in 'acc1' and 'acc2' variables. The cell structures are developed as:

    'Periods (secs)' 'RotD50 Sa (g)' 'RotD100 Sa (g)' 
