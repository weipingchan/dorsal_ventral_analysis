function wingGrid=generate4dimGrids(dimList,intv,wingGrid)
    for dimID=1:size(dimList,1)
        dim=dimList(dimID,:);
        xori=(dim(1)-1)*intv*2+1;
        yori=(dim(2)-1)*intv*2+1;

        inGrid=zeros(3,3,2);
        inGrid(1,:,:)=[wingGrid(xori,yori,:) , wingGrid(xori,yori+intv,:) , wingGrid(xori,yori+2*intv,:)];
        inGrid(2,:,:)=[wingGrid(xori+intv,yori,:) , wingGrid(xori+intv,yori+intv,:) , wingGrid(xori+intv,yori+2*intv,:)];
        inGrid(3,:,:)=[wingGrid(xori+2*intv,yori,:) , wingGrid(xori+2*intv,yori+intv,:) , wingGrid(xori+2*intv,yori+2*intv,:)];


        outGridCen=generateGrid4centroid(inGrid);
        inGridPlot1=reshape(inGrid,[],2);
        outGridPlot1=reshape(outGridCen,[],2);

        wingGrid(xori+intv/2,yori+intv/2,:)=outGridCen(1,1,:); %UL corner
        wingGrid(xori+intv/2,yori+3*intv/2,:)=outGridCen(1,2,:); %UR corner
        wingGrid(xori+3*intv/2,yori+intv/2,:)=outGridCen(2,1,:); %LL corner
        wingGrid(xori+3*intv/2,yori+3*intv/2,:)=outGridCen(2,2,:); %LR corner

        %Fill the 45 degree rotation matrix
        intvT=intv/2;
        xoriT=xori+intv;
        yoriT=yori;

        inGridT=zeros(3,3,2);
        inGridT(1,:,:)=[wingGrid(xoriT,yoriT,:)  , wingGrid(xoriT-intvT,yoriT+intvT,:) , wingGrid(xoriT-2*intvT,yoriT+2*intvT,:)];
        inGridT(2,:,:)=[wingGrid(xoriT+intvT,yoriT+intvT,:) , wingGrid(xoriT,yoriT+2*intvT,:) , wingGrid(xoriT-intvT,yoriT+3*intvT,:)];
        inGridT(3,:,:)=[wingGrid(xoriT+2*intvT,yoriT+2*intvT,:) , wingGrid(xoriT+intvT,yoriT+3*intvT,:) , wingGrid(xoriT,yoriT+4*intvT,:)];

        outGridCenT=generateGrid4centroid(inGridT);
        inGridPlotT=reshape(inGridT,[],2);
        outGridPlotT=reshape(outGridCenT,[],2);

        wingGrid(xoriT,yoriT+intvT,:)=outGridCenT(1,1,:); %Left
        wingGrid(xoriT-intvT,yoriT+2*intvT,:)=outGridCenT(1,2,:); %Up
        wingGrid(xoriT+intvT,yoriT+2*intvT,:)=outGridCenT(2,1,:); %Bottom
        wingGrid(xoriT,yoriT+3*intvT,:)=outGridCenT(2,2,:); %Right
    end
end