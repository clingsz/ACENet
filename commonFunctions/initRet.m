function r = initRet(ret,k,method)

r = ret;
Y = ret.Y;
R = ret.R;
R0 = ret.R0;

% Initialization
Phenos = size(Y,3);
Families = size(Y,2);
Twins = 2*Families;
PenPreds = size(R,3);
UnPenPreds = size(R0,3);
TotPreds = UnPenPreds + PenPreds;

ACEMask = [1e-1 1e-2 0.1];

thetas = repmat(ACEMask(:),[1 k]);
w0s = 0.0000*randn(UnPenPreds,k);
ws = 0.0001*randn(PenPreds,k);

% First Assignment for Qs

if nargin<3
    method = 'none';
end


if strcmp(method,'kmeans')
    qs = zeros(Phenos,k);
    YFlat = reshape(Y,[Twins Phenos]);
    [assn,Mus] = kmeans(YFlat',k);
    Mus = Mus';
    for i = 1:Phenos
        qs(i,assn(i)) = 1;
    end
    [~,assn] = max(qs,[],2);
else
    qs = ones(Phenos,k)+rand(Phenos,k)*0.005;
    qs = qs./repmat(sum(qs,2),[1 k]);
    Mus = randn(Twins,k)*0.01;
    [~,assn] = max(qs,[],2);
end

r.qs =qs;
r.assn =assn;
r.thetas = thetas;
r.w0s = w0s;
r.ws = ws;
r.Mus = Mus;

