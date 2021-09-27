function gridXY=deriveGridCoordinates(gridPt, tarPt)

        gridptsEdge=[reshape(gridPt(1,:,:),[],2) ; reshape(gridPt(2:end,end,:),[],2) ; reshape(gridPt(end,1:end-1,:),[],2) ; reshape(gridPt(2:end-1,1,:),[],2)];
        %gridpts=reshape(gridPt,[],2);
        close2gridPts=findClosest2Pts(gridptsEdge,tarPt);
        
        yind=round(gridPt(:,:,1)-close2gridPts(1,1),3)==0;
        xind=round(gridPt(:,:,2)-close2gridPts(1,2),3)==0;
        [gridx1, gridy1]=find(xind.*yind);
            
        yind=round(gridPt(:,:,1)-close2gridPts(2,1),3)==0;
        xind=round(gridPt(:,:,2)-close2gridPts(2,2),3)==0;
        [gridx2, gridy2]=find(xind.*yind);
        gridDist = sqrt(sum(bsxfun(@minus, close2gridPts,tarPt).^2,2));
        
        gridXY = [gridx1, gridy1]*(1-gridDist(1)/sum(gridDist))+[gridx2, gridy2]*(1-gridDist(2)/sum(gridDist)); %The coordination on grid coordinates
end