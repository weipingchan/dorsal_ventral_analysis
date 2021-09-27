function [outSeg4Pts,wingGrids ]=foreWingGrids(wingMask,dorsal_key,LRside,numberOfIntervalDegree)
% forewing=dorsal_seg==1;
% wingMask=forewing;

[B0,~]=bwboundaries(wingMask);
edgePt=B0{1};

if LRside=='L'
    lower_wing_body_join0=dorsal_key(7,:);
    keyPts0=[dorsal_key(2,:) ;dorsal_key(10,:) ];
elseif LRside=='R'
    lower_wing_body_join0=dorsal_key(9,:);
    keyPts0=[dorsal_key(3,:) ;dorsal_key(11,:) ];
end

keyPts=keyPts0;
for keyID=1:length(keyPts0)
    keyPts(keyID,:)=findCloestPt([edgePt(:,2),edgePt(:,1)],keyPts0(keyID,:)); %Added April 15, 2020
end

lower_wing_body_join=findCloestPt([edgePt(:,2),edgePt(:,1)],lower_wing_body_join0); %Added April 15, 2020

tipPtAdj=findCloestPt([edgePt(:,2),edgePt(:,1)],keyPts(2,:));
keyPts=[keyPts(1,:);tipPtAdj];% Adjust tip point to fit on the edge
shoulderMIdPt=(keyPts(1,:)+keyPts(2,:))/2;

wingAxis=diff(keyPts);
wingOrtho=reshape(null(wingAxis(:).'),1,[]);

orthoSeg=[shoulderMIdPt-wingOrtho*max(size(wingMask)) ; shoulderMIdPt+wingOrtho*max(size(wingMask))];
[intersectX,intersectY]= polyxpoly(orthoSeg(:,1),orthoSeg(:,2),edgePt(:,2),edgePt(:,1));
if length(intersectX)<2
    [intersectX,intersectY]= polyxpoly(orthoSeg(:,1)-1,orthoSeg(:,2)-1,edgePt(:,2),edgePt(:,1)); %move the line a little bit to derive two intersection points
end
wingCorner=findFarestPt([intersectX,intersectY],shoulderMIdPt);
%reg3Pts=[keyPts(1,:); keyPts(2,:); wingCorner];


%prerun to seg the wing lead into two
% kPts0=[keyPts(1,:); keyPts(2,:)];
% [~,intPts0] = edgeEvenPt(edgePt,kPts0,2);
% wingLeadMidPt=[intPts0(2,2),intPts0(2,1)];
% seg4Pts=[keyPts(1,:); wingLeadMidPt; keyPts(2,:); wingCorner];

if LRside=='L'
    inSeg4Pts=[lower_wing_body_join ; keyPts(1,:); keyPts(2,:) ; wingCorner];
elseif LRside=='R'
    inSeg4Pts=[wingCorner ; keyPts(2,:); keyPts(1,:) ; lower_wing_body_join];
end

outSeg4Pts=[lower_wing_body_join ; keyPts(1,:); keyPts(2,:) ; wingCorner];
%generate grids
%numberOfIntervalDegree=5; %minimum is 3
wingGrids=generateGrids(inSeg4Pts,edgePt,numberOfIntervalDegree);

end