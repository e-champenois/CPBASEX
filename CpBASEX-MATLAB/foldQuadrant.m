function M_out = foldQuadrant(M,x0,y0,quadrant_filter)
%
% foldQuadrant Generate a quadrant from a full image.
%
% M_out = foldQuadrant(M,x0,y0,quadrant_filter) returns the sum of the
% quadrants of the full image M specified by quadrant_filter around the
% pixel (x0,y0).
%
% Inputs:
%
% M, a 2-D array (or 3-D array for multiple images) representing the image
% to be folded.
%
% x0 and y0, integers specifying the center of the image(s).
%
% quadrant_filter, a 1-D boolean array with 4 entries specifying which
% quadrants to consider in the sum. If not specified, the default is to use
% all four quadrants.
%
% Outputs:
%
% M_out, a 2-D array (or 3-D array for multiple images) representing the
% added quadrants.

if nargin<4
    quadrant_filter = [1,1,1,1];
end

sx = size(M,2);
sy = size(M,1);
lx = min([max(sx-x0,min(~quadrant_filter([1,4]))*10^10),max(x0-1,min(~quadrant_filter([2,3]))*10^10)]);
ly = min([max(sy-y0,min(~quadrant_filter([3,4]))*10^10),max(y0-1,min(~quadrant_filter([1,2]))*10^10)]);
M_out = zeros(ly+1,lx+1,size(M,3));
quad_signs = [-1,1;-1,-1;1,-1;1,1];

for i = 1:4
    if quadrant_filter(i)
        M_out = M_out + M(y0:quad_signs(i,1):y0+quad_signs(i,1)*ly,x0:quad_signs(i,2):x0+quad_signs(i,2)*lx,:);
    end
end

end