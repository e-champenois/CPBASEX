function out = resizeFolded(data,rmax)
%
% resizeFolded Resize a quadrant folded image to a specified radius.
%
% out = resizeFolded(data,rmax) changes the maximal radius of a quadrant
% folded image by either cropping it or padding it with zeros.

if rmax>size(data,1)
    out = data;
    out(end+1:rmax,:,:)=zeros(rmax-size(data,1),size(data,2),size(data,3));
    out(:,end+1:rmax,:)=zeros(rmax,rmax-size(data,2),size(data,3));
else
    out = data(1:rmax,1:rmax,:);
end

end