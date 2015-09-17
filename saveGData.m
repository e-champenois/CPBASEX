tic
x = 0:511;
y = x;
xkratio = 8;
if xkratio == 1
    k = x;
else
    k = x(xkratio/2:xkratio:end)+0.5;
end
l = 0:2:4;
params = 0.7*xkratio;
rBF = @(x,k,params) exp(-(x-k).^2./(2*params(1)^2))./k.^2; %leave alone
zIP = @(r,k,params) sqrt((sqrt(2*10)*params(1)+k).^2-r^2); %leave alone
gData = struct('x',x,'y',y,'k',k,'l',l,'params',params,'rBF',rBF,'zIP',zIP);
G = findG(gData);
Ginv = findGinv(gData);
[U,S,V] = svd(G,0);
Up = U';
Sinv = 1./diag(S);
save(['G_r',num2str(numel(x)),'_k',num2str(numel(k)),'_l',num2str(max(l)),'_s',num2str(round(params,2)),'.mat'],'-v7.3','Up','Sinv','V','x','y','k','l','rBF','params','Ginv');
toc