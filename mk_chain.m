function [chain] = mk_chain(img)
% mk_chain Generate an 8 connected chain from an image 
%   
% input: an image path - should contain a single connected contour
% output: [Nx2] 8 connected chain of N coordinates
%

%get first connected component
[start_r,start_c] = find(img,1,'first');

%as bwtraceboundary needs an intial direction, choose the 
%first one that works
%The order make sure the outline is discovered clockwise

if img(start_r-1,start_c) == 1
    dir = 'N';
elseif img(start_r-1,start_c-1) == 1
    dir = 'NW';
elseif img(start_r-1,start_c+1) == 1
    dir = 'NE';
elseif img(start_r,start_c+1) == 1
    dir = 'E';
elseif img(start_r,start_c-1) == 1
    dir = 'W';
elseif img(start_r+1,start_c-1) == 1
    dir = 'SW';
elseif img(start_r+1,start_c+1) == 1
    dir = 'SE';
elseif img(start_r+1,start_c) == 1
    dir = 'S';
else
    assert(0);
end

chain = bwtraceboundary(img,[start_r,start_c],dir);
end