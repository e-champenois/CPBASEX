tic
x = 0:63;
xkratio = 8;
if xkratio == 1
    k = x;
else
    k = x(xkratio/2:xkratio:end)+0.5;
end
l = 0:2:4;
nk = numel(k);
nl = numel(l);
params = 0.7*xkratio;
rBF = @(x,k,params) exp(-(x-k).^2./(2*params(1)^2))./k.^2; %leave alone
zIP = @(r,k,params) sqrt((sqrt(2*10)*params(1)+k).^2-r^2); %leave alone
gData = struct('x',x,'k',k,'l',l,'params',params,'rBF',rBF,'zIP',zIP);
[K,R] = meshgrid(k,x);
frk = rBF(R,K,params);
G = findG(gData,0);
Ginv = findGinv(gData,0);
[U,S,V] = svd(G,0);
Up = U';
S = diag(S);
h5_filename = ['G_r',num2str(numel(x)),'_k',num2str(numel(k)),'_l',num2str(max(l)),'_s',num2str(round(params,2)),'.h5'];
x=x';Up=Up';S=S';V=V';frk=frk';Ginv=Ginv';
save(h5_filename,'-v7.3','x','k','l','nk','nl','Up','S','V','frk','Ginv');
toc