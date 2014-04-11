function [beta]=coordAscentENet(y,x,lambdas,alpha,init)
% check input sizes
assert(size(y,1)==size(x,1) && size(y,1)>1);
n = size(y,1);
k = size(x,2);
if length(lambdas)==1
    lambdas = lambdas*ones(k,1);
end

if ~isempty(init)
    beta = init;
else
    beta = zeros(k,1);
end

% assume default tolerance and number of iterations
TOL = 1e-5;
MAXIT = 1000;

% tracking likelihood
lls = zeros(MAXIT,1);
prevll = -realmax;


DEBUG = 0;
ll = loglik(y,x,lambdas,alpha,beta);
iter = 0;
initll = ll;
if DEBUG
    close all;
    figure(1);
    clf
    ax1 = axes('Position',[0.1 0.85 0.7 0.045])
    ax2 = axes('Position',[0.1 0.8 0.7 0.045]);
    ax3 = axes('Position',[0.81 0.3 0.09 0.49]);
    ax4 = axes('Position',[0.1 0.3 0.7 0.49]);
    ax5 = axes('Position',[0.1 0.05 0.8 0.2]);
end
x2 = zeros(k,1);
for j=1:k
    x2(j) = sum(x(:,j).^2);
end
while (ll-prevll) > TOL && iter < MAXIT
    iter = iter+1;
    prevll = ll;
    
    yres = y - x*beta;
    
    % updates    
    for j=1:k
        yres = yres + beta(j)*x(:,j);
        beta(j) = 1/(x2(j) + lambdas(j)*(alpha)) * shrinkThreshold(yres'*x(:,j),(1-alpha)*lambdas(j));
        yres = yres - beta(j)*x(:,j);
    end
    
    % likelihood for new state
    ll = loglik(y,x,lambdas,alpha,beta);
    
    assert(ll-prevll>=-1e-10)
    
    if DEBUG
        lls(iter) = ll;
        
        
        axes(ax1);cla;
        imagema(y',0);
        title(['ElasticNet fit with $$\sum_i (y_i - \beta_0 - x_i''\beta)^2 - '...
            num2str(alpha*mean(lambdas(2:end))/2) '\sum_j \beta_j^2 - ' ...
            num2str((1-alpha)*mean(lambdas(2:end))) '\sum_j |\beta_j|$$'],'Interpreter','Latex');
        
        set(gca,'XTick',[],'YTick',1,'YTickLabel','y','TickLength',[0 0]);
        
        axes(ax2);cla;
        imagema((x*beta)',0);
        set(gca,'XTick',[],'YTick',1,'YTickLabel','X*beta','TickLength',[0 0]);
        
        axes(ax3);cla;
        barh(flipud(beta));
        title('beta')
        set(gca,'YTick',[]);
        ylim([0.5 length(beta)+0.5])
        xlim([-max(abs(beta)) max(abs(beta))])
        
        axes(ax4);cla;
        imagema(x',0);
        set(gca,'XTick',[],'YTick',[size(x,2)/2],'YTickLabel','X');
        
        axes(ax5);cla;
        plot(lls(1:iter))
        xlabel('iteration');
        ylabel('log-likelihood');
        
        drawnow
    end
end
if iter == MAXIT
    fprintf('coordAscentENet Algorithm did not converge: %d\n',ll-prevll);
%     error('converge')
end

function ll = loglik(y,x,lambdas,alpha,beta)
n = size(y,1);
ll = -1/2*(y - x*beta)'*(y - x*beta)...
    -alpha/2*sum(lambdas.*beta.^2) - (1-alpha)*sum(abs(lambdas.*beta));

function b = shrinkThreshold(b,lambda)
b = sign(b).*max(abs(b) - lambda,0);

