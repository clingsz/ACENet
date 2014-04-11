function logProb = comPheno(y,w,R,invL2,Mu,ALL,mp)

if nargin<5
    Mu = mean(y,2);
end

if nargin<6
    ALL = 0;
end

if nargin<7
    mp = ones(size(y,2),1);
end

logProb = zeros(size(y,2),1);

zyg = ones(length(Mu)/2,1);

sigma = 0.3;

invL1 = computeCholeskyInvCov([0 0 sigma],zyg);
r1 = Mu - R*w;
x = invL1*r1;
val1 = sum(log(diag(invL1)))-1/2*((x')*x);
%val1 = log p(Mu|w0,sigma)

for i=1:size(y,2)
    % val2 = log p(y_i|Mu,theta)
    r2 = y(:,i) - Mu;
    x = invL2*r2;
    val2 = sum(log(diag(invL2)))-1/2*((x')*x);
%     disp(size(mp));
    logProb(i) = val2 + log(mp(i));
end

if ALL==1
    logProb = sum(logProb) + val1;
end
