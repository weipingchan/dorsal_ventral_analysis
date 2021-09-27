function segLine_dorsal=projectSegLine(wingMask,in_key,segLine,LRside)
    if LRside=='L'
        keyPts=[in_key(7,:) ; in_key(1,:)];
    elseif LRside=='R'
        keyPts=[in_key(9,:) ; in_key(4,:)];
    end
%     
%     if nnz(wingMask)==0 %If one wing is missing, provide a fake one
%         fakepolygon=[in_key(3,:); keyPts; in_key(11,:); ];
%         wingMask=roipoly(wingMask,fakepolygon(:,1),fakepolygon(:,2));
%     end
    
    [B0,~]=bwboundaries(wingMask);
    edgePt=flip(B0{1},2);

    wing_body_join=findCloestPt([edgePt(:,1),edgePt(:,2)],keyPts(1,:));
    fore_hind_wing_join=findCloestPt([edgePt(:,1),edgePt(:,2)],keyPts(2,:));
        
    [tformSeg,~,~] = estimateGeometricTransform([segLine(end,:) ; segLine(1,:)] , [wing_body_join ; fore_hind_wing_join], 'similarity');
    [X0,Y0] = transformPointsForward(tformSeg,segLine(:,1),segLine(:,2));
    segLine_dorsal=[X0,Y0];

end