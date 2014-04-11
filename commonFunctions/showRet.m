function showRet(ret,FL)

% maxMod = ret.para.maxMod;

if nargin<2
    FL = 0;
end

KK=0;
if KK==1
pthetas = ret.thetas(:,ret.assn);
ret.pws = ret.ws(:,ret.assn);
act = sum(abs(ret.pws)>1e-1)';

[~,nlst] = sort(pthetas(1,:));

tit = {'pThetas','activeSet'};
out = {pthetas(:,nlst),act(nlst)};

pl = length(tit);
figure(99+FL);
for i =1:pl
    if size(out{i},2)>1
        subplot(pl,1,i),imagesc(out{i}),title(tit{i}),colorbar;
    else
        subplot(pl,1,i),bar(out{i}),title(tit{i});
    end
end
end
%% Show Module
SHOWMODULE = 1;
if SHOWMODULE
    pN = 4;
    
    k = size(ret.ws,2);
    Phenos = size(ret.Y,3);
    thetas = ret.thetas;
    w0s = ret.w0s;
    ws = ret.ws;
    Y = ret.Y;
    
    
    if ~isfield(ret,'qs')
        assn = ret.assn;
        qs = zeros(Phenos,k);
        for i = 1:Phenos
            qs(i,assn(i))=1;
        end
    else
        qs =ret.qs;
    end
    
    lst = find(sum(qs,2));
    [~,assn] = max(qs,[],2);
    bucket = hist(assn,k);
    showthetas = thetas;
    showthetas(:,bucket==0)=0;
    
    KK = 1;
    
    YFlat = reshape(Y,[2*size(Y,2) size(Y,3)]);
    figure(100+FL);
    subplot(pN/KK,KK,1),imagesc(qs),ylabel('qs'),colorbar;
%     subplot(pN/KK,KK,2),bar(showthetas','stack'),ylabel('thetas');
    subplot(pN/KK,KK,2),imagesc(thetas),ylabel('thetas'),colorbar;
%     subplot(pN/KK,KK,3),bar(sum(abs(ws)>0)),ylabel('Ws');
    subplot(pN/KK,KK,3),imagesc(ws),ylabel('ws'),colorbar;
%     subplot(pN/KK,KK,4),imagesc(YFlat),ylabel('Y'),colorbar;
    TT = YFlat;
    PP = [];
    
    for i = 1:size(qs,2)
        PP = [PP TT(:,assn==i)];
    end
    subplot(pN/KK,KK,4),imagesc(PP),ylabel('Y sort'),colorbar;
%     subplot(pN/KK,KK,6),hist(assn,k),ylabel('Bucket'),colorbar;
    
end
% figure(101+FL),imagesc(PP),colorbar;