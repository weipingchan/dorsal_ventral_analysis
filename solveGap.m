function wingGriddeep=solveGap(wingGriddeep)
    %solve the final empty points
    [x0,y0]=find(wingGriddeep(:,:,1)==0);
    holelocList=[x0,y0];

    for holeID=1:size(holelocList,1)
        holeloc=holelocList(holeID,:);
        cendat=[wingGriddeep(holeloc(1)-1,holeloc(2),:) ; wingGriddeep(holeloc(1),holeloc(2)+1,:) ; wingGriddeep(holeloc(1)+1,holeloc(2),:) ; wingGriddeep(holeloc(1),holeloc(2)-1,:)];
        cen=mean(cendat);
        inPt = inpolygon(cen(1),cen(2),cendat(:,:,1),cendat(:,:,2)); %Test if the centoid is in the polygon or not
        if inPt==0 %If not, use the closest point as the centroid
            cen=findCloestPt(reshape(cendat,[],2),reshape(cen,[],2));
        end
        wingGriddeep(holeloc(1),holeloc(2),:)=cen;
    end

end