function tailinfo=tail_module2(gridsParameter_dorsal, scale, tailMinArea, side)
    if side=='L'
        refineAreaH=gridsParameter_dorsal{5}{1};
        originalAreaH=gridsParameter_dorsal{7};
        gridPt=gridsParameter_dorsal{3}{2};
        seg4Pts=gridsParameter_dorsal{3}{1};
    elseif side=='R'
        refineAreaH=gridsParameter_dorsal{6}{1};
        originalAreaH=gridsParameter_dorsal{8};
        gridPt=gridsParameter_dorsal{4}{2};
        seg4Pts=gridsParameter_dorsal{4}{1};
    end
%figure,imshowpair(originalAreaH,refineAreaH);

    tailraw=immultiply(originalAreaH,imcomplement(refineAreaH));
    tailParts=bwareaopen(imdilate(imerode(tailraw,strel('disk',2)),strel('disk',2)),tailMinArea); %The value here determine how small can a tail be defined
    forkRef=seg4Pts(3,:);
    bodyRef=seg4Pts(1,:);

    [tailL,nParts]=bwlabel(tailParts);

    disp(['Totaly find ',num2str(nParts),' tail regions.']);
    minBranchLength=5;
    if nParts>0 %run if we have at least one tail
        %BWtail0 = bwskel(logical(originalAreaH),'MinBranchLength',minBranchLength);
        BWtail0 = bwskel(logical(originalAreaH),'MinBranchLength',minBranchLength);
        
        branchimage0 = bwmorph(BWtail0, 'branchpoints');
        ss = regionprops(branchimage0,'centroid');
        forkPts0=cat(1,ss.Centroid);

        try
            forkPt=findCloestPt(forkPts0,forkRef);
        catch
            [ip,jp] = find(bwmorph(immultiply(BWtail0,imcomplement(tailParts)),'endpoints'));
            forkPt0=findCloestPt([ip,jp],bodyRef);
            forkPt=flip(forkPt0); %modified Sep, 30, 2020
            disp('No fork point. Use the endpoint close to body as the fork point.');
        end

        tailinfo=cell(0,nParts);
        for tailID=1:nParts
            try
                tail=tailL==tailID;

                [i,j] = find(immultiply(bwmorph(BWtail0,'endpoints'),tail));
                tailStat=regionprops(tail,'Eccentricity');

                if ~isempty(i) && tailStat.Eccentricity>0.75
                    Dd = bwdistgeodesic(BWtail0,forkPt(1),forkPt(2),'quasi');

                    %Measure the length from ref point to the top of tail
                    distanceTail=[];
                    for n = 1:numel(i)
                        distanceTail=[distanceTail; [i(n),j(n),Dd(i(n),j(n))]];
                    end
                    tailTip2Fork0 = distanceTail(distanceTail(:,3)==max(distanceTail(:,3)),:);
                    tailTip2Fork=tailTip2Fork0(3);
                    tailTip=tailTip2Fork0(:,1:2);

                    %Measure the length from ref point to the base of tail         
                    [i2,j2] = find(immultiply(bwmorph(immultiply(BWtail0,imcomplement(tailParts)),'endpoints'),imdilate(tail,strel('disk',3))));
                    distanceTail2=[i2,j2,Dd(i2,j2)];
                    tailBase2Fork0 = distanceTail2;
                    tailBase2Fork = tailBase2Fork0(3);
                    tailBase=tailBase2Fork0(1:2);

%                     figure,imshowpair(tail,BWtail0)
%                     hold on;
%                     plot(forkPt(1),forkPt(2),'rO')
%                     plot(tailTip(2),tailTip(1),'rX')
%                     plot(tailBase(2),tailBase(1),'rX')

                    tailElongation=tailTip2Fork-tailBase2Fork+minBranchLength;
                    tailArea=nnz(tail);
                    tailWavg=tailArea/tailElongation;
                    tailCurvature=(tailElongation-minBranchLength)/pdist([tailBase;tailTip ]);
                    
                else
                    boundary0=immultiply(imdilate(tail,strel('disk',1)),refineAreaH);
                    boundaryIt0=bwskel(boundary0);
                    [bi, bj]=find(boundaryIt0);
                    intPts0 = interparc(3,bi,bj);
                    tailBase=intPts0(2,:);

                    boundaryAreaB = bwboundaries(imdilate(boundaryIt0,strel('disk',3))); 
                    boundaryBpts= boundaryAreaB{1};
                    tailB = bwboundaries(tail); 
                    tailBpts= tailB{1};
                    inBoundaryArea = inpolygon(tailBpts(:,1),tailBpts(:,2),boundaryBpts(:,1),boundaryBpts(:,2));
                    tarPts=tailBpts(~inBoundaryArea,:);
                    tailDist = sqrt(sum(bsxfun(@minus, tarPts,tailBase).^2,2));
                    tailElongation=max(tailDist);
                    tailTip=tarPts(tailDist==max(tailDist),:);
                    tailTip=tailTip(1,:);

                    tailArea=nnz(tail);
                    tailWavg=tailArea/tailElongation;
                    tailCurvature=tailElongation/sqrt(sum(bsxfun(@minus, tailTip, tailBase).^2,2));
                end
                    %Measure the length of the shared boundary between tail and
                    %main wing part
                    boundary=immultiply(imdilate(tail,strel('disk',1)),imdilate(refineAreaH,strel('disk',1)));
                    boundaryIt=bwskel(boundary);
        %             stat=regionprops(boundaryIt,'perimeter');  %faster way
        %             tailBoundaryL = (stat.Perimeter-2)/2;  %34.4220
                    [ib,jb] = find(bwmorph(boundaryIt,'endpoints'));
                    DdB = bwdistgeodesic(boundaryIt,jb(1),ib(1),'quasi');
                    tailBoundaryL = DdB(ib(2),jb(2));

                %%
                %Find the corresponding grid for a given tail Base and boundary
                %length
                tailBasegridXY=deriveGridCoordinates(gridPt, flip(tailBase));
                tailBaseEndgridXY1=deriveGridCoordinates(gridPt, flip([ib(1),jb(1)]));
                tailBaseEndgridXY2=deriveGridCoordinates(gridPt, flip([ib(2),jb(2)]));

                tailBaseRaw=[[ib(1),jb(1)]; tailBase; [ib(2),jb(2)]];
                tailBaseGrid=[tailBaseEndgridXY1; tailBasegridXY; tailBaseEndgridXY2];

                tailMeasurement=[tailBoundaryL/scale, tailElongation/scale, tailArea/scale^2, tailWavg/scale, tailCurvature];

                tailinfo{tailID}= {tail, tailBaseRaw, tailBaseGrid, tailMeasurement};
                %{}= tail part
                %{}{1}=tail mask
                %{}{2}=XY coordination on the original image
                %{}{3}=XY coordination projected on the summarized grid (begin from 1)
                %{}{4}= [tail boundary Length, tail length, tail area, tail width, tailCurvature; curvature is unit less
            catch
                disp(['No. ',num2str(tailID),' tail failed in processing.']);
            end
        end
        disp(['All ',num2str(nParts),' tail regions have been successfully processed.']);
    else
        tailinfo={-9999};    
    end

end