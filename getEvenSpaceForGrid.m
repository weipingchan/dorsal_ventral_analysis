function evenIntEdge=getEvenSpaceForGrid(seg4Pts,edgePt,numberOfIntervalDegree)
    numberOfInterval=2^numberOfIntervalDegree;

    kPts=seg4Pts(1:2,:);
    [segEdgePts1,intPts1] = edgeEvenPt2(edgePt,kPts,numberOfInterval);
    ointPts1=flip([intPts1(1,:);flip(intPts1(2:length(intPts1)-1,:),1);intPts1(length(intPts1),:)],2); %reoder points
    %1-> begining point
    %2-> in order (from end point to begining)
    %last ->end point

    kPts=[seg4Pts(3,:); seg4Pts(2,:)];
    [segEdgePts2,intPts2] = edgeEvenPt2(edgePt,kPts,numberOfInterval);
    ointPts2=flip([intPts2(1,:);intPts2(2:length(intPts2)-1,:);intPts2(length(intPts2),:)],2); %reoder points

    kPts=[seg4Pts(4,:); seg4Pts(3,:)];
    [segEdgePts3,intPts3] = edgeEvenPt2(edgePt,kPts,numberOfInterval);
    ointPts3=flip([intPts3(1,:);intPts3(2:length(intPts3)-1,:);intPts3(length(intPts3),:)],2); %reoder points

    kPts=[seg4Pts(4,:); seg4Pts(1,:)];
    [segEdgePts4,intPts4] = edgeEvenPt2(edgePt,kPts,numberOfInterval);
    ointPts4=flip([intPts4(1,:);flip(intPts4(2:length(intPts4)-1,:),1);intPts4(length(intPts4),:)],2); %reoder points

    evenIntEdge={ointPts1, ointPts2, ointPts3, ointPts4};
end