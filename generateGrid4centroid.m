function outGridCen=generateGrid4centroid(inGrid)
% inGrid=zeros(3,3,2);
% inGrid(1,:,:)=[wingGrid(1,1,:);wingGrid(1,numberOfInterval/2+1,:) ; wingGrid(1,numberOfInterval+1,:)];
% inGrid(2,:,:)=[wingGrid(numberOfInterval/2+1,1,:);wingGrid(numberOfInterval/2+1,numberOfInterval/2+1,:) ; wingGrid(numberOfInterval+1,numberOfInterval+1,:)];
% inGrid(3,:,:)=[wingGrid(numberOfInterval+1,1,:);wingGrid(numberOfInterval+1,numberOfInterval/2+1,:) ; wingGrid(numberOfInterval+1,numberOfInterval+1,:)];

outGridCen=zeros(2,2,2);

% itr=2;
%UL corner
cendat=[inGrid(1,1,:) ; inGrid(1,2,:) ; inGrid(2,2,:) ; inGrid(2,1,:)];
cen=mean(cendat);
inPt = inpolygon(cen(1),cen(2),cendat(:,:,1),cendat(:,:,2)); %Test if the centoid is in the polygon or not
if inPt==0 %If not, use the closest point as the centroid
    cen=findCloestPt(reshape(cendat,[],2),reshape(cen,[],2));
end
outGridCen(1,1,:)=cen;

%UR corner
cendat=[inGrid(1,3,:) ; inGrid(2,3,:) ; inGrid(2,2,:) ; inGrid(1,2,:)];
cen=mean(cendat);
inPt = inpolygon(cen(1),cen(2),cendat(:,:,1),cendat(:,:,2)); %Test if the centoid is in the polygon or not
if inPt==0 %If not, use the closest point as the centroid
    cen=findCloestPt(reshape(cendat,[],2),reshape(cen,[],2));
end
outGridCen(1,2,:)=cen;

%LL corner
cendat=[inGrid(3,1,:) ; inGrid(2,1,:) ; inGrid(2,2,:) ; inGrid(3,2,:)];
cen=mean(cendat);
inPt = inpolygon(cen(1),cen(2),cendat(:,:,1),cendat(:,:,2)); %Test if the centoid is in the polygon or not
if inPt==0 %If not, use the closest point as the centroid
    cen=findCloestPt(reshape(cendat,[],2),reshape(cen,[],2));
end
outGridCen(2,1,:)=cen;

%LR corner
cendat=[inGrid(3,3,:) ; inGrid(3,2,:) ; inGrid(2,2,:) ; inGrid(2,3,:)];
cen=mean(cendat);
inPt = inpolygon(cen(1),cen(2),cendat(:,:,1),cendat(:,:,2)); %Test if the centoid is in the polygon or not
if inPt==0 %If not, use the closest point as the centroid
    cen=findCloestPt(reshape(cendat,[],2),reshape(cen,[],2));
end
outGridCen(2,2,:)=cen;

end