function ret = stanTwinData(Y)
% Standardize Twin Data
% input matrix Y ~ [2 F P]

STANDARIZE = 1;

    F = size(Y,2);
    P = size(Y,3);
    
    y = reshape(Y,[F*2 P]);
    
    for i = 1:P
        x = y(:,i);
        x = (x-mean(x));
        x = x/norm(x);
        if STANDARIZE==1
            x = x*sqrt(F*2-1);
        end
        y(:,i) = x;
    end
       
    ret = reshape(y,[2 F P]);
end