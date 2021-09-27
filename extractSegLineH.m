function  segLine=extractSegLineH(wingMask,ventral_key_flip,LRside)
    %Use left wing to develop the code
%     wingMask=wingLH;
%     LRside='L';
    if LRside=='L'
        keyPts=[ventral_key_flip(7,:) ; ventral_key_flip(1,:)];
    elseif LRside=='R'
        keyPts=[ventral_key_flip(9,:) ; ventral_key_flip(4,:)];
    end

    [B0,~]=bwboundaries(wingMask);
    edgePt=flip(B0{1},2);

    wing_body_join=findCloestPt([edgePt(:,1),edgePt(:,2)],keyPts(1,:));
    fore_hind_wing_join=findCloestPt([edgePt(:,1),edgePt(:,2)],keyPts(2,:));

    iniPt=find(edgePt(:,1)==wing_body_join(1) & edgePt(:,2)==wing_body_join(2));
    if LRside=='R'
        loopEdgePt=[edgePt(iniPt:size(edgePt,1),:) ; edgePt(1:iniPt-1,:)];
        endPt=find(loopEdgePt(:,1)==fore_hind_wing_join(1) & loopEdgePt(:,2)==fore_hind_wing_join(2));
        segLine=flip(loopEdgePt(1:endPt,:),1);
    elseif LRside=='L'
        loopEdgePt=[edgePt(iniPt+1:size(edgePt,1),:) ; edgePt(1:iniPt,:)];
        endPt=find(loopEdgePt(:,1)==fore_hind_wing_join(1) & loopEdgePt(:,2)==fore_hind_wing_join(2));
        segLine=loopEdgePt(endPt:length(loopEdgePt),:);
    end

end