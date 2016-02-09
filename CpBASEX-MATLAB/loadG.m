function gData = loadG(filename, full)

%info = h5info(filename);
%dsets = info.Datasets;
% for i = 1:numel(dsets)
%     name = dsets(i).Name;
%     gData.(name) = h5read(filename,['/',name]);
% end
if nargin<2 || full
    gData.Ginv = h5read(filename,'/Ginv');
end
gData.x = h5read(filename,'/x');
gData.y = h5read(filename,'/y');
gData.k = h5read(filename,'/k');
gData.l = h5read(filename,'/l');
gData.params = h5read(filename,'/params');
gData.rBF = @(x,k,params) exp(-(x-k).^2./(2*params(1)^2))./k.^2;
gData.Up = h5read(filename,'/Up');
gData.Sinv = h5read(filename,'/Sinv');
gData.V = h5read(filename,'/V');

end