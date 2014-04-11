function r = getModuleMean(ret)

k = size(ret.ws,2);
thetaMean = zeros(3,k);
assn = ret.assn;
Y = ret.Y;
zyg = ret.zyg;
ws = ret.ws;
thetas = ret.thetas;

parfor m = 1:k
    disp(m);
    if (sum(assn==m)>0)
       theta=findThetaMean(Y(:,:,assn==m),zyg);
       if sum(ret.ws(:,m))==0
           theta2 = [0 0 0]';
       else
            Pred = getPred(ret,m);
            Pred = reshape(Pred,[2 size(Y,2) size(Pred,2)]);
            theta2 = findThetaMean(Pred,zyg);
       end
    end
    disp('module count activeSet thetas');
    disp([m sum(assn==m) sum(abs(ws(:,m))>0) theta']);
    tm{m} = theta;
    tr{m} = theta2;
end
for m = 1:k
    if (sum(assn==m)>0)
        thetaMean(:,m)=tm{m};
        thetaRW(:,m)=tr{m};
    end
end
    

pthetas = thetas(:,assn);
pws = ws(:,assn);

r = ret;
r.thetaMean = thetaMean;
r.thetaRW = thetaRW;
r.pthetas = pthetas;
r.pws = pws;
