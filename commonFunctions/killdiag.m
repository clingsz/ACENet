function P = killdiag(P)
N = size(P,1);
for i = 1:N
    P(i,i)=0;
end
