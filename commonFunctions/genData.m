function [x zyg] = genData(theta,Families,repl,seed)
% example: [x zyg] = genData([1 1 1],50,2,1);
% x: 2 X N


if nargin>=4
    randn('seed',seed);
end

if nargin<3
    repl = 1;
end




a = theta(1); c = theta(2); e = theta(3);

%Families = 200;
N = Families*repl;
%N = Families;

zyg = repmat([ones(Families/2,1); -ones(Families/2,1);],[repl,1]);
%zyg = repmat([-ones(Families/2,1); ones(Families/2,1);],[1,1]);

ind = 0.5*(zyg + 1) + 1; % zyg=1 -> ind=2, zyg=-1 -> ind=1
[Smz,Sdz] = computeCov([a,c,e]);
Ss(:,:,1) = Sdz;
Ss(:,:,2) = Smz;
x = zeros(2,N);
for i=1:N
    S = Ss(:,:,ind(i)); % S is the cov matrix for the ith twin
    L = chol(S);  % L is the cholvsky decomposition of S
    r = randn(1,2); % r is the random number pair, mean is 0
    vec = (r*L);
    x(1:2,i) = vec;    
end

%x = repmat(x,[1,repl]);
%zyg = repmat(zyg,[repl,1]);

% function test
% [x,z] = genData(theta,T);
% [Smz,Sdz] = computeCov(theta);
% cov(x(:,find(z==1)))-Smz;
% cov(x(:,find(z==-1)))-Sdz;
% end
