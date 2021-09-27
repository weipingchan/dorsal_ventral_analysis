function [segEdgePts,intPts] = edgeEvenPt2(edgePt,kPts,Nsec)
kPt1=flip(kPts(1,:));%Change coordination for numeric calculation
kPt2=flip(kPts(2,:));%Change coordination for numeric calculation

closePts=[findClosest2Pts(edgePt,kPt1); findClosest2Pts(edgePt,kPt2)];
subDist=pdist2(closePts(1:2,:), closePts(3:4,:));
subLoc=find(subDist==min(min(subDist)));
if subLoc==1
    innerPts=[closePts(1,:); closePts(3,:)];
elseif subLoc==2
    innerPts=[closePts(2,:); closePts(3,:)];
elseif subLoc==3
    innerPts=[closePts(1,:); closePts(4,:)];
else
    innerPts=[closePts(2,:); closePts(4,:)];
end

[~, rLoc1] = ismember(innerPts(1,:),edgePt, 'rows');
[~, rLoc2] = ismember(innerPts(2,:),edgePt, 'rows');

segEdgePts0=edgePt(min([rLoc1,rLoc2]):max([rLoc1,rLoc2]),:);
segEdgePts1=[edgePt(max([rLoc1,rLoc2]):end,:) ; edgePt(1:min([rLoc1,rLoc2]),:)];

%Calculate the distance of two different segment to see which one is closer
%to our key points
Dist01 = sum(sqrt(sum(bsxfun(@minus, segEdgePts0, kPt1).^2,2)));
Dist02 = sum(sqrt(sum(bsxfun(@minus, segEdgePts0, kPt2).^2,2)));
Dist11 = sum(sqrt(sum(bsxfun(@minus, segEdgePts1, kPt1).^2,2)));
Dist12 = sum(sqrt(sum(bsxfun(@minus, segEdgePts1, kPt2).^2,2)));
Dist0=Dist01+Dist02;
Dist1=Dist11+Dist12;

if Dist0>Dist1
    segEdgePts2=segEdgePts1;
else
    segEdgePts2=segEdgePts0;
end

if pdist2(kPt1,segEdgePts2(1,:))<pdist2(kPt2,segEdgePts2(1,:))
    segEdgePts3=[kPt1; segEdgePts2; kPt2] ;
elseif pdist2(kPt1,segEdgePts2(1,:))>pdist2(kPt2,segEdgePts2(1,:))
    segEdgePts3=[kPt2; segEdgePts2; kPt1] ;
end

%Remove duplicates but keep the order
[~,AdjInd] = unique(segEdgePts3,'rows'); 
segEdgePts4=segEdgePts3(sort(AdjInd),:);

%Randomly remove n points from the segment in order to make 'interparc'
%work. (added April 18, 2020)
deleterown=0;
while 1
    firstRow=segEdgePts4(1,:);
    endRow=segEdgePts4(end,:);
    shrinkData=segEdgePts4(2:end-2,:);
    deleteidx=randsample(1:size(shrinkData,1),deleterown) ;
    shrinkData(deleteidx,:)=[]; % remove those rows
    segEdgePts=[firstRow; shrinkData; endRow];
    try
        intPts0 = interparc(Nsec+1,segEdgePts(:,1),segEdgePts(:,2));
        break
    catch
        deleterown=deleterown+1;
    end
    if deleterown>size(shrinkData,1)/4
        break
    end
end

%Nsec=16;
% intPts0 = interparc(Nsec+1,segEdgePts(:,1),segEdgePts(:,2));
intPts=[kPt1 ;intPts0(2:end-1,:); kPt2];
end