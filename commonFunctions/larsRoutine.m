function [beta,activeSet] = larsRoutine(X,y,method,DEBUG)
% Find the beta using lars and significance test
% X ~ N by P
% Y ~ N by 1
% beta ~ P by 1
% activeSet ~ P by 1

if nargin<3
    method = 'sig';
end

if nargin<4
DEBUG = 0;
end
% sig for significant test
% bic for bic test

% X = X * 10;
% y = y * 10;

[N,P] = size(X);
assert(N==size(y,1));


% y = y-mean(y);
% y = y/norm(y)*sqrt(N-1);

lambda = realmax;
activeSet = [];
alpha = expinv(0.95);
lambdas = realmax*ones(P,1);

maxDOF = 10;

betas = larsen(X,y);

if strcmp(method,'sig')
    warning off;
%     keyboard;
    CR = 10;
    if (N<=P)
        az = abs(betas(CR,:))>0;
        Xin = X(:,az);
        sigma2 = norm(y-Xin*(Xin\y))^2/(N-CR);
    else
        sigma2 = norm(y-X*(X\y))^2/(N-P);
    end
    
    bestDOF =1;
    ffk = zeros(1,maxDOF);
    
    for DOF=2:maxDOF
        beta = betas(DOF,:)';
        lambdas(activeSet) = lambda;
        [betaA]=coordAscentENet(y,X,lambdas,0,[]);
        F1 = y'*(X*beta);
        F2 = y'*(X*betaA);
        Fk = F1 - F2;
        Fk = Fk/sigma2;

        ff1(DOF)=F1;
        ff2(DOF)=F2;
        ffk(DOF)=Fk;
%         keyboard;
        if Fk>alpha
            %             fprintf('Fk small!!  %.2f  < %.2f\n',Fk,alpha);
            bestDOF = DOF;
        else
            break;
        end
        lambda = max(abs(X'*(y-X*beta)));
        activeSet = find(beta);
        
    end
    
    if DEBUG
        alst = 0.95:0.005:0.995;
        for j = 1:length(alst)
            aa(j,:) = expinv(alst(j))*ones(1,maxDOF);
        end
        subplot(2,1,1),plot([ffk;aa]'),legend('fk');
        subplot(2,1,2),plot(betas(1:maxDOF,:));
%         keyboard;
    end
%     [~,bestDOF] = max(ffk);
    beta = betas(bestDOF,:)';
elseif strcmp(method,'bic')
%     scale = sqrt(N-1);
%     y = y-mean(y);
%     y = y/norm(y)*scale;
%     for i = 1:P
%         xx = X(:,i);
%         xx = xx-mean(xx);
%         xx = xx./norm(xx)*scale;
%         X(:,i)=xx;
%     end
    
    bic = realmax*ones(1,maxDOF);
    for DOF=1:maxDOF
        beta = betas(DOF,:)';
        r = y-X*beta;
        ll(DOF) = -5*norm(r);
        bic(DOF) = aicbic(ll(DOF),DOF,N);
    end
    [~,bestDOF] = min(bic);
    beta = betas(bestDOF,:)';
    DOF = bestDOF;
    activeSet = find(beta);
    dl = [1:maxDOF]*log(N);
    if DEBUG
        subplot(2,1,1),plot([bic;-2*ll;dl]'),legend('bic','-ll','dl');
        subplot(2,1,2),plot(betas(1:maxDOF,:));
    end
    beta = betas(DOF,:)';
end


