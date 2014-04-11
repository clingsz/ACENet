function in = createSyn(Modules)
% This is a function to create a synthetic input data
% output variable in will be inputable to fitHM.m
%
% INPUT:
%
% Modules: the modules that want to be created
%
% OUTPUT:
%
% in is a structure contains following attributes:
%
% Y : 2 x F x Phenos, the twin phenotype measurements
% R : 2 x F x PenPreds, the penalized twin Regulator measurements
% R0 : 2 x F x UnPenPreds, the unpenalized twin regulator measuremetns
% zyg: F x 1, the zygosity for each familiy (1 for monozygotic, -1 for dizygotic twins)
%

disp('Creating Synthetic Data...');
ActiveSet = 5;

R0reds = 1;
Preds = 20+R0reds;

Families = 60;
Phenos = 400;
Twins = 2*Families;

Rreds = Preds-R0reds;

Y = zeros(2,Families,Phenos);
R = ones(1,Families*Preds);


RR = R(Families*R0reds+1:end);
R0 = R(1:Families*R0reds);

R = [R;R];
RR = [RR;RR];
R0 = [R0;R0];

R = reshape(R,[2,Families,Preds]);
RR = reshape(RR,[2,Families,Rreds]);
R0 = reshape(R0,[2,Families,R0reds]);


w = zeros(1,Modules);
w0s = repmat(w,[R0reds,1]);
ws = repmat(w,[Rreds,1]);

for m = 1:Modules
    Rlst = randperm(Rreds,ActiveSet);
    ws(Rlst,m) = 0.5.*linspace(1,0.1,ActiveSet).*sign(randn(1,ActiveSet));
end


thetas = [linspace(0,2,Modules)',0.1*ones(Modules,1),1*ones(Modules,1)];
thetas = 0.2*sqrt(thetas.^2./repmat(sum(thetas.^2,2),[1,3]));


for m = 1:Modules
    Pm = Rreds/Modules;
    [x,zyg] = genData([0 0 1],Families,Pm);  % c[2*Families] z[Families*1]
    eps = reshape(x,[Twins Pm]);
    RRFlat(:,(m-1)*Pm+1:Pm*m)=eps;
end
%ws = [zeros(1,Modules);ws];

%w = 0.1*randn(Preds,Modules);

YFlat = reshape(Y,[Twins Phenos]);
R = reshape(RRFlat,[2 Families Rreds]);
RFlat = [reshape(R0,[Twins R0reds]) reshape(R,[Twins Rreds])];
% keyboard;
w = [w0s;ws];

mu = RFlat*w; % [Twins*Modules]
qs = zeros(Phenos,Modules);

for m = 1:Modules
    Pm = Phenos/Modules;
    assn((m-1)*Pm+1:Pm*m) = m;
    [x,zyg] = genData(thetas(m,:),Families,Pm);  % c[2*Families] z[Families*1]
    eps = reshape(x,[Twins Pm]);
    c = repmat(mu(:,m),[1,Pm])+eps;
    YFlat(:,(m-1)*Pm+1:Pm*m) = c;
    qs((m-1)*Pm+1:Pm*m,m) = 1;
end

zyg = zyg(1:Families)';
Y = reshape(YFlat,[2,Families,Phenos]);

ref.Modules = Modules;
ref.thetas = thetas';
ref.w0s = w0s;
ref.ws = ws;
ref.assn = assn;
ref.Mus = mu;
ref.Y = Y;
ref.R = R;
ref.R0 = R0;
ref.zyg = zyg;

in.Y = Y;
in.R = R;
in.R0 = R0;
in.zyg = zyg(1,:);
in.ref = ref;

