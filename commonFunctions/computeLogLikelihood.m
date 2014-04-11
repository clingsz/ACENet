function [ll,dtheta,ddtheta] = computeLogLikelihood(x,mzyg,intheta,q,fast)
USEOLD=0;
if USEOLD
    [ll,dtheta,ddtheta] = computelld(x,mzyg,intheta);
else
    oldx = x;
    %keyboard;
    
    if nargin<4
        q = ones(2,size(x,2));
    end
    
    
    if nargin>=5
        x = fast.x;
        nmz = fast.nmz;
        ndz = fast.ndz;
        xmz = fast.xmz;
        xdz = fast.xdz;
        xmzxmzT = fast.xmzxmzT;
        xdzxdzT = fast.xdzxdzT;
        
    else
        x = x.*sqrt(q);
        nmz = sum(q(1,mzyg==1));
        ndz = sum(q(1,mzyg==-1));
        xmz = x(:,mzyg==1);
        xdz = x(:,mzyg==-1);
        xmzxmzT = xmz*xmz';
        xdzxdzT = xdz*xdz';
    end
    
    theta = zeros(3,1); dtheta = zeros(3,1); ddtheta = zeros(3,1);
    theta(:) = intheta(:);
    O2 = ones(2,2);
    I2 = eye(2);
    Amz = O2;                Cmz = O2; Emz = I2;
    Adz = (1/2*O2 + 1/2*I2); Cdz = O2; Edz = I2;
    
    a = theta(1); c = theta(2); e = theta(3);
    
    Xmz(:,:,1) = O2; Xmz(:,:,2) = O2; Xmz(:,:,3) = I2;
    Xdz(:,:,1) = (1/2*O2 + 1/2*I2); Xdz(:,:,2) = O2; Xdz(:,:,3) = I2;
    
    
    Smz = zeros(2,2); Sdz = zeros(2,2);
    for i=1:3
        Smz = Smz + theta(i).^2*Xmz(:,:,i);
        Sdz = Sdz + theta(i).^2*Xdz(:,:,i);
    end
    
    oSmz = a^2*Amz + c^2*Cmz + e^2*Emz;
    oSdz = a^2*Adz + c^2*Cdz + e^2*Edz;
    
    assert(max(max(abs(oSmz - Smz)))<1e-10);
    assert(max(max(abs(oSdz - Sdz)))<1e-10);
    % assert(max(eigs(oSmz))>-1e-10);
    % assert(max(eigs(oSdz))>-1e-10);
    
    
    iSmz = inv(Smz);
    iSdz = inv(Sdz);
    
    %nmz = sum(mzyg==1);
    %ndz = sum(mzyg==-1);
    %     keyboard;
    
    
    
    ll = - nmz*log(2*pi) - nmz/2*log(det(Smz)) - 1/2*sum(sum((xmzxmzT).*iSmz)) ...
        - ndz*log(2*pi) - ndz/2*log(det(Sdz)) - 1/2*sum(sum((xdzxdzT).*iSdz));
    
    
    
    Mmz = -nmz*iSmz + iSmz'*xmzxmzT*iSmz';
    Mdz = -ndz*iSdz + iSdz'*xdzxdzT*iSdz';
    odtheta = zeros(3,1);
    odtheta(1) = a*sum(sum(Mmz.*Amz + Mdz.*Adz));
    odtheta(2) = c*sum(sum(Mmz.*Cmz + Mdz.*Cdz));
    odtheta(3) = e*sum(sum(Mmz.*Emz + Mdz.*Edz));
    
    
    for i=1:3
        dtheta(i) = theta(i)*sum(sum(Mmz.*Xmz(:,:,i) + Mdz.*Xdz(:,:,i)));
    end
    
    if ~(max(abs(dtheta - odtheta))<1e-10)
        disp('calculate dtheta error');
        disp(dtheta);
%         keyboard;
    end
    
    
    
    Namz = nmz*iSmz*(2*a*Amz)*iSmz - 2*iSmz*(2*a*Amz)*iSmz'*xmzxmzT*iSmz';
    Ncmz = nmz*iSmz*(2*c*Cmz)*iSmz - 2*iSmz*(2*c*Cmz)*iSmz'*xmzxmzT*iSmz';
    Nemz = nmz*iSmz*(2*e*Emz)*iSmz - 2*iSmz*(2*e*Emz)*iSmz'*xmzxmzT*iSmz';
    
    
    Nadz = ndz*iSdz*(2*a*Adz)*iSdz - 2*iSdz*(2*a*Adz)*iSdz'*xdzxdzT*iSdz';
    Ncdz = ndz*iSdz*(2*c*Cdz)*iSdz - 2*iSdz*(2*c*Cdz)*iSdz'*xdzxdzT*iSdz';
    Nedz = ndz*iSdz*(2*e*Edz)*iSdz - 2*iSdz*(2*e*Edz)*iSdz'*xdzxdzT*iSdz';
    
    for i=1:3
        M2mz(:,:,i) = nmz*iSmz*(2*theta(i)*Xmz(:,:,i))*iSmz - 2*iSmz*(2*theta(i)*Xmz(:,:,i))*iSmz'*xmzxmzT*iSmz';
        M2dz(:,:,i) = ndz*iSdz*(2*theta(i)*Xdz(:,:,i))*iSdz - 2*iSdz*(2*theta(i)*Xdz(:,:,i))*iSdz'*xdzxdzT*iSdz';
    end
    
    oddtheta = zeros(3,1);
    oddtheta(1) = sum(sum(Mmz.*Amz + Mdz.*Adz)) + a*sum(sum(Namz.*Amz + Nadz.*Adz));
    oddtheta(2) = sum(sum(Mmz.*Cmz + Mdz.*Cdz)) + c*sum(sum(Ncmz.*Cmz + Ncdz.*Cdz));
    oddtheta(3) = sum(sum(Mmz.*Emz + Mdz.*Edz)) + e*sum(sum(Nemz.*Emz + Nedz.*Edz));
    
    for i=1:3
        for j=1:3
            ddtheta(i,j) = theta(i)*sum(sum(M2mz(:,:,j).*Xmz(:,:,i) + M2dz(:,:,j).*Xdz(:,:,i)));
            if i==j
                ddtheta(i,j) = ddtheta(i,j) + sum(sum(Mmz.*Xmz(:,:,i) + Mdz.*Xdz(:,:,i)));
            end
        end
    end
    
    
    if ~(all(abs(diag(ddtheta) - oddtheta)<1e-10))
        keyboard;
    end
    
    % dtheta = dtheta/size(x,2);
    % ddtheta = ddtheta/size(x,2);
    
    DEBUGFORKEANS = 0;
    
    if DEBUGFORKEANS == 1
        [l1,d1,dd1] = computelld(oldx,mzyg,intheta);
        dif1 = l1-ll;
        dif2 = d1-dtheta;
        dif3 = dd1-ddtheta;
        if abs(dif1)>0
            %keyboard;
        end
    end
end