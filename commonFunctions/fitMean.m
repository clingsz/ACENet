function [w,w0,theta,invL,newMu]=fitMean(Y,R,R0,zyg,Mu)

REFIT = 1;
sigma = 0.3;

Phenos = size(Y,3);
Families = size(Y,2);
Twins = 2*Families;
PenPreds = size(R,3);
UnPenPreds = size(R0,3);
TotPreds = PenPreds + UnPenPreds;
Families = Twins/2;

if (any([size(Y,1) size(Y,2) size(Y,3)] ~= [2 Families Phenos])), error('Input Y should be of size 2 x Families x Phenos');end;
if (any([size(R,1) size(R,2) size(R,3)] ~= [2 Families PenPreds])),error('Input R should be of size 2 x Families x PenPredictors');end;
if (any([size(R0,1) size(R0,2) size(R0,3)] ~= [2 Families UnPenPreds])),error('Input R0 should be of size 2 x Families x PenPredictors');end;
if (length(zyg) ~= Families), error('Input zyg should be of length Families'); end;

familywise = @(x) reshape(x,[2 length(x)/2]);
repeatForAllPhenos = @(v) repmat(v,[Phenos 1]);
unroll = @(v) v(:);

YFlatT = reshape(Y,[2*Families Phenos])'; % Phenos * Twins
YFlat = YFlatT';
YMean = mean(YFlat,2);

yall = unroll(YFlat);

R = reshape(R,[2*Families PenPreds])';
R0 = reshape(R0,[2*Families UnPenPreds])';

zygT = repeatForAllPhenos(zyg)'; % Families x Phenos
zygAll = unroll(zygT);


w = zeros(UnPenPreds + PenPreds,1);
RBoth = [R0;R];
RFlat = RBoth';
oldll = -realmax;
newll = -realmax/2;

theta = [0 0 0.1]';

it = 0;
maxit = 20;
while (it<maxit && newll-oldll>1e-2)
    it = it +1;
    muall = repmat(Mu,[Phenos 1]);
    allPhenosFamilywiseMeanResidual = familywise(yall-muall);
    
    proj = zeros(PenPreds,UnPenPreds);
    invL = computeCholeskyInvCov([0 0 sigma],zyg);
    g = Mu'*invL';
    S = RBoth*invL';
    S0 = S(1:UnPenPreds,:);
    [res,beta0] = findRes(S0',g');
    yhat = res;
    SS = S(UnPenPreds+1:end,:);
    Xhat = zeros(2*Families,PenPreds);
    for i=1:PenPreds
        [resi,proj(i,:)]=findRes(S0',SS(i,:)');
        Xhat(:,i) = resi;
    end
    
    invL = computeCholeskyInvCov(theta,zyg);
    prevW = w;
    beforeWAll = comPheno(YFlat,w,RFlat,invL,Mu,1);
    
    
    beta = larsRoutine(Xhat,yhat,'sig');
    
    w = [beta0;beta];
    
    if REFIT ==1
        strongLambdas = 1000*ones(length(w),1);
        strongLambdas(abs(w)>1e-6) = 0;
        [w] = coordAscentENet(g',S',strongLambdas,0,w);
    end
    
    afterWAll = comPheno(YFlat,w,RFlat,invL,Mu,1);
    
    change = afterWAll-beforeWAll;
    if change<-1e-3
        %fprintf('broken W update');
        w = prevW;
        afterWAll = beforeWAll;
        disp(change);
    end
    
    
    
    if Phenos > 1
        theta = newton(allPhenosFamilywiseMeanResidual,zygAll,[1 1 1]);
    else
        theta = [0 0 1]';
    end
    
    
    invL = computeCholeskyInvCov(theta,zyg);
    afterTheta = comPheno(YFlat,w,RFlat,invL,Mu,1);
    change = afterTheta - afterWAll;
    if change<-1e-3
        fprintf('broken Theta update');
        disp(change);
    end
    
    invS = invL'*invL;
    invS0 = 1/sigma^2*diag(ones(Twins,1));
    mu0 = RBoth'*w;
    
    
    newMu = (invS0 + Phenos*invS)\(Phenos*invS*YMean+invS0*mu0);
    
    
    
    afterMu = comPheno(YFlat,w,RFlat,invL,newMu,1);
    change = afterMu - afterTheta;
    if change<-1e-3
        fprintf('broken Mu update');
        disp(change);
        newMu = Mu;
    end
    Mu = newMu;
    oldll = newll;
    newll = afterMu;
end

if (it==maxit)
    disp('reach maxit!');
end

w0 = w(1:UnPenPreds);
w = w(UnPenPreds+1:end);

