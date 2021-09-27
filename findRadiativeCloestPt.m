function cloestPt=findRadiativeCloestPt(ptList,tarPt)
    cen=mean(ptList);
    vec = [cen(1)-tarPt(1),cen(2)-tarPt(2)];
    segLine = [tarPt+2*vec; tarPt-2*vec];
    
    [xi,yi] = polyxpoly(ptList(:,1),ptList(:,2),segLine(:,1),segLine(:,2));
    
    intersectPts=[xi,yi];

    if size(intersectPts,1)==1
        cloestPt=intersectPts;
    else
        cloestPt=findCloestPt(intersectPts,tarPt);
    end       
end