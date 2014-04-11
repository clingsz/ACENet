function Pred = getPred(ret,m)

R0 = reshape(ret.R0,[2*size(ret.R0,2) size(ret.R0,3)]);
R = reshape(ret.R,[2*size(ret.R,2) size(ret.R,3)]);
RFlat = [R0 R];
w = [ret.w0s(:,m);ret.ws(:,m)];
Pred = RFlat*w;