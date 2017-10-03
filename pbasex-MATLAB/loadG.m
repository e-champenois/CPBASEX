function gData = loadG(filename, make_images)

gData.x = h5read(filename,'/x')';
gData.nk = h5read(filename,'/nk');
gData.nl = h5read(filename,'/nl');
gData.Up = h5read(filename,'/Up')';
gData.S = h5read(filename,'/S')';
gData.V = h5read(filename,'/V')';
gData.frk = h5read(filename,'/frk')';

if nargin==2 && make_images
    gData.Ginv = h5read(filename,'/Ginv')';
end

end