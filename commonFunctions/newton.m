function [theta,success] = newton(x,zyg,init,q,DEBUG)
% x 2 by Family by 1
% zyg 2 by Family
% q Pheno by 1

if nargin<5
    DEBUG = 0;
end

if nargin<4
    q = ones(size(x));
end

if nargin<3
    init = [1 1 1];
end


theta = 0.1*ones(3,1);
if exist('init','var')
    theta(:) = init;
end

ntheta = zeros(3,1);

usable = find(init);
[ll,d,dd] = computeLogLikelihood(x,zyg,theta,q);
it = 0;
c1 = 1e-4; c2 = 0.9;
[e,v] = eigs(dd(usable,usable));

if any(diag(v)>0)
    dd(usable,usable) = e*diag(min(diag(v),-1e-2))*e';
end

p = ones(length(usable),1);
accepted = 1;

fast.x = x.*sqrt(q);
fast.xmz = fast.x(:,zyg==1);
fast.nmz = sum(q(1,zyg==1));
fast.xdz = fast.x(:,zyg==-1);
fast.ndz = sum(q(1,zyg==-1));
fast.xmzxmzT = fast.xmz*fast.xmz';
fast.xdzxdzT = fast.xdz*fast.xdz';

while norm(p)>1e-4 && it<5000 && norm(d(usable))>1e-4
    accepted = 0;
    step = 1;
    while ~accepted && step>1e-20        
        p = (dd(usable,usable))\(-d(usable));        
        
        ntheta(usable) = theta(usable) + step*p;
        
        
        [nll,nd,ndd] = computeLogLikelihood(x,zyg,ntheta,q,fast);
        if nll <= ll + c1*step*p'*d(usable) || ...
           (abs(p'*d(usable)) > 1e-10 && abs(p'*nd(usable)) <= c2*abs(p'*d(usable)))
            step = step*0.75;
        else
            accepted = 1;
        end
    end
    [e,v] = eigs(ndd(usable,usable));
    if any(diag(v)>0)
        ndd(usable,usable) = e*diag(min(diag(v),-1e-2))*e';
    end
    
    d = nd;
    dd = ndd;
    ll = nll;
    theta = ntheta;
    it = it+1;
    if DEBUG       
        lls(it) = ll;
        theta_hist(it,:) = theta;
        if (mod(it,100)==0)
        subplot(2,1,1),plot(lls);
        subplot(2,1,2),plot(theta_hist),legend('a','c','e');
        
        drawnow
        disp('it ll normd');
        disp([it ll norm(d(usable))]);
        end
    end
    
    if ~accepted
%         keyboard;
        break;
    end
end
if ~accepted
    fprintf('Newton algorithm did not converge: %d!\n',norm(d(usable)));
    %error('convergence')
end
success = accepted;

theta = abs(theta);
