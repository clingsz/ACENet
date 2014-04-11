function Yres = getRes(Y,R,w)

% Y 2 F P
% R 2 F R
% w R 1

Phenos = size(Y,3);
Families = size(Y,2);
Twins = 2*Families;
PenPreds = size(R,3);
Families = Twins/2;

YFlat = reshape(Y,[Twins Phenos]);
RFlat = reshape(R,[Twins PenPreds]);

Yres = YFlat - repmat(RFlat*w,[1 Phenos]);
Yres = reshape(Yres,[2 Families Phenos]);

