function  segLine=extractSegLineF(wingMask,dorsal_key,LRside)
    %Use left wing to develop the code
    %wingMask=wingLF;
    if LRside=='L'
        keyPts=[dorsal_key(7,:) ; dorsal_key(1,:)];
    elseif LRside=='R'
        keyPts=[dorsal_key(9,:) ; dorsal_key(4,:)];
    end
% 
%     if nnz(wingMask)==0 %If one wing is missing, provide a fake one
%         fakepolygon=[dorsal_key(3,:); keyPts; dorsal_key(11,:); ];
%         wingMask=roipoly(wingMask,fakepolygon(:,1),fakepolygon(:,2));
%     end
        
    [B0,~]=bwboundaries(wingMask);
    edgePt=flip(B0{1},2);

    wing_body_join=findCloestPt([edgePt(:,1),edgePt(:,2)],keyPts(1,:));
    fore_hind_wing_join=findCloestPt([edgePt(:,1),edgePt(:,2)],keyPts(2,:));

    iniPt=find(edgePt(:,1)==wing_body_join(1) & edgePt(:,2)==wing_body_join(2));
    if LRside=='L'
        loopEdgePt=[edgePt(iniPt:size(edgePt,1),:) ; edgePt(1:iniPt-1,:)];
        endPt=find(loopEdgePt(:,1)==fore_hind_wing_join(1) & loopEdgePt(:,2)==fore_hind_wing_join(2));
        segLine=flip(loopEdgePt(1:endPt,:),1);
    elseif LRside=='R'
        loopEdgePt=[edgePt(iniPt+1:size(edgePt,1),:) ; edgePt(1:iniPt,:)];
        endPt=find(loopEdgePt(:,1)==fore_hind_wing_join(1) & loopEdgePt(:,2)==fore_hind_wing_join(2));
        segLine=loopEdgePt(endPt:length(loopEdgePt),:);
    end

end