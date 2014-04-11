function [b,LL,DOF] = getBIC(r)
% need: Y R0 R w0s ws zyg assn thetas

bound = 0;

DOF = sum(sum(abs(r.ws)>bound)) + sum(sum(abs(r.thetas)>bound))+sum(sum(abs(r.Mus)>bound));

Y = r.Y;
R0 = r.R0;
R = r.R;
w0s = r.w0s;
ws = r.ws;
zyg = r.zyg;
assn = r.assn;
thetas = r.thetas;
Mus = r.Mus;
YFlat = reshape(Y,[size(Y,2)*2 size(Y,3)]);
RFlat = reshape(cat(3,R0,R),[size(R,2)*2 size(R0,3)+size(R,3)]);
mps = r.mps;

LL = 0;

maxMod = max(assn);

for m = 1:maxMod
    if sum(assn==m)>0
        theta = thetas(:,m);
        w0 = w0s(:,m);
        w = ws(:,m);
        Mu = Mus(:,m);
        mp = mps(:,m);
        invL = computeCholeskyInvCov(theta,zyg);
        lnow = comPheno(YFlat(:,assn==m),[w0;w],RFlat,invL,Mu,1,mp);
        LL = LL + sum(lnow);
    end
end

[~,b] = aicbic(LL,DOF,size(Y,2)*2);
