function [im] = texture_checker_board(m,n,s,color)

x = 1:m;
y = 1:n;
[X,Y] = meshgrid(x,y);

% s = 16;

F = mod(ceil((X)/s) + ceil((Y)/s),2);
% 0.3, 0.3, 0.8

im = ones(m,n,3);
im(:,:,1) = color(1) * 255 * F;
im(:,:,2) = color(2) * 255 * F;
im(:,:,3) = color(3) * 255 * F;

im = uint8(im);

