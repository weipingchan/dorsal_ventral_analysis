function [outSeg4Pts,wingGrids ]=hindWingGrids(wingMask,inKey,inRegKey,LRside,numberOfIntervalDegree)
% inKey=ventral_key_flip;
% inRegKey=regPtLH;
% wingMask=refineAreaLH;

[B0,~]=bwboundaries(wingMask);
edgePt=B0{1};

wing_body_join=inRegKey(1,:);
if LRside=='L'
    keyPts0=[inKey(6,:) ;inRegKey(2,:) ];
elseif LRside=='R'
    keyPts0=[inKey(5,:) ;inRegKey(2,:) ];
end

keyPts=keyPts0;
for keyID=1:length(keyPts0)
    keyPts(keyID,:)=findCloestPt([edgePt(:,2),edgePt(:,1)],keyPts0(keyID,:)); %Added April 15, 2020
end

tipPtAdj=findCloestPt([edgePt(:,2),edgePt(:,1)],keyPts(2,:));
keyPts=[keyPts(1,:);tipPtAdj];% Adjust tip point to fit on the edge
wingCorner=inRegKey(3,:);
%reg3Pts=[keyPts(1,:); keyPts(2,:); wingCorner];


%prerun to seg the wing lead into two
% kPts0=[keyPts(1,:); keyPts(2,:)];
% [~,intPts0] = edgeEvenPt(edgePt,kPts0,2);
% wingLeadMidPt=[intPts0(2,2),intPts0(2,1)];
% seg4Pts=[keyPts(1,:); wingLeadMidPt; keyPts(2,:); wingCorner];

if LRside=='L'
    inSeg4Pts=[keyPts(1,:) ; wing_body_join ; wingCorner ; keyPts(2,:)];
elseif LRside=='R'
    inSeg4Pts=[keyPts(2,:) ; wingCorner ; wing_body_join ; keyPts(1,:)];
end

outSeg4Pts=[wing_body_join ; keyPts(1,:); keyPts(2,:) ; wingCorner];
%generate grids
%numberOfIntervalDegree=5; %minimum is 3
wingGrids=generateGrids(inSeg4Pts,edgePt,numberOfIntervalDegree);

end