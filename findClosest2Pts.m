function closePts=findClosest2Pts(edgePt,kPt)
%kPt1=reg3Pts(1,:);
edgePt2 = unique(edgePt, 'rows'); %Remove duplicates
Dist = sqrt(sum(bsxfun(@minus, edgePt2, kPt).^2,2));
sortDist=sort(Dist);
closePts=[edgePt2(Dist==sortDist(1),:);edgePt2(Dist==sortDist(2),:)];
if size(closePts,1)>2
    closePts=closePts(1:2,:);
end
end