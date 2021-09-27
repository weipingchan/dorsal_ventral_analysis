function [x,y]=orderPtsPoly(x,y)
    %Step 1: Find the unweighted mean of the vertices:
    cx = mean(x);
    cy = mean(y);
    
    %Step 2: Find the angles:
    a = atan2(y - cy, x - cx);
    
    %Step 3: Find the correct sorted order:
    [~, order] = sort(a);
    
    %Step 4: Reorder the coordinates:
    x = x(order);
    y = y(order);
end