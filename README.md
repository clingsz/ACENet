ACENet
=======

Fit the gene and regulator expressions of twins into modules

File descriptions:

example.m:
A synthetic example is shown in example.m. It generate a synthetic data set with 5 modules of phenotypes. And it will fit the phenotypes into 5 modules and estimate their regulator weights, mean and ACE variance.

addpaths.m
This script will automatically add the required folders into the current work space.

Inside the main folder:

createSyn.m
This function can generate synthetic data set with given module number.

fitHM.m
This function can fit the phenotypes and regulators in to different modules, and estimate their regulator activities, phenotype expression mean of each module, and the ACE parameter for each module.

reFit.m
This is the mean function to fit the model. Switch of Mix portions updating is in this file.
