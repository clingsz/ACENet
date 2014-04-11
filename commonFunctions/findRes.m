function [r,beta] = findRes(X,y)
% Find the residual for best regression
% r = y - beta*X;
% y ~ N * 1
% X ~ N * P
% beta ~ P * 1
% keyboard;
    beta = X\y;
    res = y - X*beta;
%     test = res*X';
%     if (sum(test>1e-5)>0)
%         disp('beta is not orthogonal');
%         disp(res);
%         save('resDbg.mat');
%     end
    r = res;
end

    