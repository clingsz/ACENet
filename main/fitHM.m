function ret = fitHM(in,maxMod)
% Fit a phenotype and regulator measurements into different modules and
% estimate their mean and variance using mixture model.
% 
% INPUT:
% in is a structure contains following attributes:
% Y : 2 x F x Phenos, the twin phenotype measurements
% R : 2 x F x PenPreds, the penalized twin Regulator measurements
% R0 : 2 x F x UnPenPreds, the unpenalized twin regulator measuremetns
% zyg: 2 x F, the zygosity for each familiy (1 for monozygotic, -1 for dizygotic twins)
%
% maxMod is the max module count allowed in the mixture model
%
% OUTPUT:
%
% ret is a structure contains following attributes:
%
% Y : 2 x F x Phenos, the twin phenotype measurements
% R : 2 x F x PenPreds, the penalized twin Regulator measurements
% R0 : 2 x F x UnPenPreds, the unpenalized twin regulator measuremetns
% zyg: F x 1, the zygosity for each familiy (1 for monozygotic, -1 for dizygotic twins)
%
% ws: PenPreds x Modules: the weights for each penalized regulators of each
% modules
% w0s: UnPenPreds x Modules: the weight for each unpenalized regulators of
% each modules
% thetas: 3 x Modules: the residual ACE parameter
% assn: Phenos x 1: The assigned module number for each phenotype
% Mus: 2F x Modules: the estimated mean for each twin in each module 
%
% Author: Tianxiang Gao
% Email: tgao@cs.unc.edu
% Release: 1.0
% Release Date: 10/27/12

% Copy the inputs
r.Y = in.Y;
r.R = in.R;
r.R0 = in.R0;
r.zyg = in.zyg;

% initialize the model
disp('Initializing Model...');
r = initRet(r,maxMod);

% fit the model
ret = reFit(r,30);





