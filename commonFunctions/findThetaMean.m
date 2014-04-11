function [theta]=findThetaMean(Y,zyg)

Families = size(Y,2);
Phenos = size(Y,3);

ACEMask = [1 1 1];

familywise = @(x) reshape(x,[2 length(x)/2]);
repeatForAllPhenos = @(v) repmat(v,[Phenos 1]);
unroll = @(v) v(:);
zygT = repeatForAllPhenos(zyg)'; % Families x Phenos
zygAll = unroll(zygT);

yFlat = reshape(Y,[2*Families Phenos]); % Twins x Phenos
% SPECIAL - twin mean
% yFlat = yFlat - repmat(mean(yFlat,2),[1 Phenos]);

yall = unroll(yFlat);
ymean= mean(mean(yFlat));
yall = yall - ymean;
yRes = familywise(yall);


theta = ACEMask;

if (size(Y,3)>0)
    qr = ones(2,Families*Phenos);
    theta = newton(yRes,zygAll,theta,qr,0);
else
    q = ones(2,Families);
    y = Y-mean(mean(Y));
    theta = newton(y,zyg,theta,q,0);
end


