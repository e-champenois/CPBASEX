%%
x = 0:10;%0:1:511;
y = x;
k = 0.5:10;%0.5:2:510.5;
l = 0:2:4;
rBF = @(x,k,params) exp(-(x-k).^2./(2*params(1)^2))/k^2;
params = 1;
%%
tic
G = findG(x,y,k,l,rBF,params);
toc
%%
tic
Ginv = findGinv(x,y,k,l,rBF,params);
toc
%%
tic
[U,S,V] = svd(G,0);
toc
Up = U';
S = diag(S);
%%
tic
save(['PBASEX_Cart/G_r',num2str(numel(x)),'_k',num2str(numel(k)),'_l',num2str(max(l)),'_s',num2str(params),'.mat'],'-v7.3','Ginv','Up','S','V','x','y','k','l','rBF','params');
toc