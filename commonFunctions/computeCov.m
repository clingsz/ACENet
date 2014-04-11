function [Smz,Sdz] = computeCov(theta)
% computeCov computes a covariance matrix specified by ACE parameters.
%
% [Smz,Sdz] = computeCov(theta)
%     theta  - 3 x 1 array containing values a,c,e
% Returns
%     Smz    - covariance for a monozygotic twin pair
%     Sdz    - covariance for a dizygotic twin pair
%
% Smz and Sdz can be expressed using standard matrices for capturing
% additive (A), common/family environmental effects (C), and independent
% environmental effects (E)
%     Smz = a^2*Amz + c^2*C + d^2*E
%     Sdz = a^2*Adz + c^2*C + d^2*E
% where we note that A matrix is dependent on zygosity of the twin pair.
%
O2 = ones(2,2);
I2 = eye(2);
Amz = O2;                Cmz = O2; Emz = I2;
Adz = (1/2*O2 + 1/2*I2); Cdz = O2; Edz = I2;

a = theta(1); c = theta(2); e = theta(3);

Smz = a^2*Amz + c^2*Cmz + e^2*Emz;
Sdz = a^2*Adz + c^2*Cdz + e^2*Edz;
