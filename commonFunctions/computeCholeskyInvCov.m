function [invL,L,Sigma] = computeCholeskyInvCov(theta,zyg)
[Smz,Sdz] = computeCov(theta);
Sigma = kron(spdiag(zyg==1),Smz) + kron(spdiag(zyg==-1),Sdz);
L = chol(Sigma,'lower');
invL = inv(L);