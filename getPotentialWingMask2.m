function [potentialWingMask,nullArea]=getPotentialWingMask2(wingMask,segLine,part)
     [B0,~]=bwboundaries(wingMask);
    edgePt=flip(B0{1},2);

    wing_body_join=findCloestPt([edgePt(:,1),edgePt(:,2)],segLine(1,:));
    fore_hind_wing_join=findCloestPt([edgePt(:,1),edgePt(:,2)],segLine(end,:));

    iniPt=find(edgePt(:,1)==wing_body_join(1) & edgePt(:,2)==wing_body_join(2));
    if part=='RH' | part=='RF'
        loopEdgePt=[edgePt(iniPt:size(edgePt,1),:) ; edgePt(1:iniPt-1,:)];
        endPt=find(loopEdgePt(:,1)==fore_hind_wing_join(1) & loopEdgePt(:,2)==fore_hind_wing_join(2));
        segLine2=flip(loopEdgePt(1:endPt,:),1);
    elseif  part=='LH' | part=='LF'
        loopEdgePt=[edgePt(iniPt+1:size(edgePt,1),:) ; edgePt(1:iniPt,:)];
        endPt=find(loopEdgePt(:,1)==fore_hind_wing_join(1) & loopEdgePt(:,2)==fore_hind_wing_join(2));
        segLine2=loopEdgePt(endPt:length(loopEdgePt),:);
    end

    segJoin=[segLine; segLine2];
    potArea=roipoly(wingMask,segJoin(:,1),segJoin(:,2));
    potMask0=imfill(wingMask+potArea,'hole');

    if part=='LH'
        rmAreaPt=[[1,1] ; segLine ; [size(wingMask,2), round(size(wingMask,1)/2)] ; [size(wingMask,2),1]];
    elseif part=='RH'
        rmAreaPt=[[size(wingMask,2),1] ; segLine ; [1,round(size(wingMask,1)/2)] ;  [1,1]];
    elseif part=='LF'
        rmAreaPt=[[1, size(wingMask,2)] ; segLine ; [size(wingMask,2), round(size(wingMask,1)/2)] ; flip(size(wingMask))];
    elseif part=='RF'
        rmAreaPt=[flip(size(wingMask)) ; segLine ; [1,round(size(wingMask,1)/2)] ; [1, size(wingMask,2)] ];
    end
    rmArea=roipoly(wingMask,rmAreaPt(:,1),rmAreaPt(:,2));
    potentialWingMask=imdilate(imerode(bwareafilt(immultiply(potMask0,imcomplement(rmArea))>0,1),strel('disk',1)),strel('disk',1)); %Remove unnecessary segment lines due to bias
    nullArea=potentialWingMask~=(potentialWingMask & wingMask);
end