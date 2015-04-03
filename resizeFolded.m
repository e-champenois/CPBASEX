function out = resizeFolded(data,rmax)
%
% resizeFolded Resize a quadrant folded image to a specified radius.
%
% out = resizeFolded(data,rmax) changes the maximal radius of a quadrant
% folded image through cropping or padding with zeros.

if size(data,1)>rmax
    x1 = rmax;
    x2 = 0;
else
    x1 = size(data,1);
    x2 = rmax-x1;
end
if size(data,2)>rmax
    y1 = rmax;
    y2 = 0;
else
    y1 = size(data,2);
    y2 = rmax-y1;
end
z = size(data,3);

out = [data(1:x1,1:y1,:),zeros(x1,y2,z);zeros(x2,rmax)];

end