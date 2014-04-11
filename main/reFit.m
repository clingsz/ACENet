function r = reFit(r,maxit)
% fit the current model with more iterations until the model converges
%
% INPUT:
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
% maxit: is the max iteration allowed for the fitting
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
% SWITCH
% MIXP = 1 calcualte with mixing portions
%
% Author: Tianxiang Gao
% Email: tgao@cs.unc.edu
% Release: 1.1
% Release Date: 10/27/12

DEBUG = 0;
MIXP = 0;

k = size(r.qs,2);
Y = r.Y;
R = r.R;
R0 = r.R0;
zyg = r.zyg;
qs = r.qs;

Phenos = size(Y,3);
Twins = size(Y,2)*2;
YFlat = reshape(Y,[2*size(Y,2) Phenos]);
RFlat = [reshape(R0,[2*size(R0,2) size(R0,3)]) reshape(R,[2*size(R,2) size(R,3)])];

Mus = r.Mus;
ws = r.ws;
w0s = r.w0s;
thetas = r.thetas;
logProbs = zeros(size(Y,3),k);

if nargin<2
    maxit = 1;
end

mps = ones(Phenos,k)*1/k;

[~,assn] = max(qs,[],2);
qs = zeros(Phenos,k);
for i = 1:Phenos
    qs(i,assn(i))=1;
end
if MIXP == 1
    for m = 1:k
        mpval = sum(assn==k)/Phenos;
        if mpval<1e-6
            mpval = 1/Phenos;
        end
        mps(:,k) = mpval;
    end
end
r.mps = mps;


oldll = -realmax;
[~,nowll,~] = getBIC(r);

disp('Start Model Fitting...');

for it = 1:maxit
    if abs(oldll-nowll)<1e-4
        break;
    end
    
    for i =1:k
        parY{i} = Y(:,:,assn==i);
        if (sum(assn==i)<=2)
            buc = hist(assn,1:k);
            [~,maxBuc] = max(buc);
            randPhenoLst = find(assn==maxBuc);
            randPhenoLst = randPhenoLst(randperm(length(randPhenoLst)));
            randPheno = randPhenoLst(1:3);
            parY{i} = Y(:,:,randPheno);
            Mus(:,i) =  zeros(Twins,1);
        end
    end
    
    for i=1:k
        disp([num2str(it) ' iteration: Fitting for module ' num2str(i)]);
        [w,w0,theta,invL,Mu]=fitMean(parY{i},R,R0,zyg,Mus(:,i));
        parws{i} = w;
        parw0s{i} = w0;
        parthetas{i} = theta;
        parMus{i} = Mu;
        parlogProbs{i} = comPheno(YFlat,[w0;w],RFlat,invL,Mu,0,mps(:,i));
    end
    
    for i=1:k
        ws(:,i) = parws{i};
        w0s(:,i) = parw0s{i};
        thetas(:,i) = parthetas{i};
        logProbs(:,i) = parlogProbs{i};
        Mus(:,i) = parMus{i};
    end
    
    logProbs = logProbs - repmat(max(logProbs,[],2),[1 k]);
    qs = exp(logProbs);
    qs = qs./repmat(sum(qs,2),[1 k]);
    
    [~,assn] = max(qs,[],2);
    
    qs = zeros(Phenos,k);
    for i = 1:Phenos
        qs(i,assn(i))=1;
    end
    
    if MIXP == 1
        for m = 1:k
            mpval = sum(assn==k)/Phenos;
            if mpval<1e-6
                mpval = 1/Phenos;
            end
            mps(:,k) = mpval;
        end
    end
    
    r.mps = mps;
    r.w0s =w0s;
    r.ws = ws;
    r.thetas = thetas;
    r.assn =assn;
    r.qs=qs;
    r.Mus = Mus;
    
    oldll = nowll;
    [bic,nowll,~] = getBIC(r);
    
    fprintf('%.0f Iteration get a model with BIC = %f\n',it,bic);
    %     disp(bic);
    if DEBUG
        disp(nowll-oldll);
        save('data/temp.mat','r');
        showRet(r);
        drawnow;
    end
    
end
if (it==maxit)
    disp('maxit reached!');
end
disp('Fitting finished!');
