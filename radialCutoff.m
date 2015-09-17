function mat = radialCutoff(mat,x0,y0,radius)

[X,Y] = meshgrid(1:size(mat,2),1:size(mat,1));
R = sqrt(abs(X-x0).^2+abs(Y-y0).^2);
mat(repmat(R,1,1,size(mat,3))>radius) = 0;

end