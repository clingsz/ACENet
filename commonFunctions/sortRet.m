function nret = sortRet(ret,nlst)

% Things related to module number:
% assn,thetaMean,ws,w0s,thetas,qs
nret = ret;
ret.thetaMean = ret.thetas;
if nargin<2
    [~,nlst] = sort(ret.thetaMean(1,:),'descend');
end
nret.ws = ret.ws(:,nlst);
nret.w0s = ret.w0s(:,nlst);
nret.thetas = ret.thetas(:,nlst);
nret.thetaMean = ret.thetaMean(:,nlst);
nret.qs = ret.qs(:,nlst);
nret.Mus = ret.Mus(:,nlst);

nret.assn = ret.assn;
for i =1:length(ret.assn)
    nret.assn(i) = find(nlst==ret.assn(i));
end