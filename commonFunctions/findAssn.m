function assn = findAssn(R,ret)

k = size(ret.ws,2);

logProb = zeros(k,1);
zyg = ones(size(R,2),1);

r = reshape(R,[2*size(R,2) size(R,3)]);

for i = 1:k
    theta = ret.thetas(:,i);
    Mu = ret.Mus(:,i);
    invL2 = computeCholeskyInvCov(theta,zyg);
    % val2 = log p(y_i|Mu,theta)
    r2 = r - Mu;
%     keyboard;
    x = invL2*r2;
    val2 = sum(log(diag(invL2)))-1/2*((x')*x);
    logProb(i) = val2;
end

[~,assn] = max(logProb);

    