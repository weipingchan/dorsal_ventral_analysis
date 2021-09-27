function [refineArea,regPtH,reconstructRegPtH]=refineHindWing2(wingMask,in_key,LRside,ornamentRatio)
    %Use wingLH as model for development
%     wingMask=wingLH;
%     [B0,~]=bwboundaries(wingMask);
%     edgePt=B0{1};

    %Remove ornament
    amountErode=round(nnz(wingMask)/3000);
    wingreduce=imdilate(bwareafilt(imerode(wingMask,strel('disk',amountErode)),1),strel('disk',amountErode));

    %make a fully closed shape
    chain = mk_chain(flip(wingreduce)); %Flip the mask to counter the flipping effect in the following process

    [chainCode] = chaincode(chain); %Derive the chain code of the shape
    chainBeginPt=[chainCode.x0,chainCode.y0];
    maskChainCode=transpose(chainCode.code);

    % Calculte traversl distance
    xx = calc_traversal_dist(maskChainCode);
    % Starting point is assumed from [0 0]
    reconstShape0 = [0 0; xx]+chainBeginPt;
    reconstShape=[reconstShape0(:,1),size(wingreduce,1)-reconstShape0(:,2)];

    basicHarmonic=5; %Derive basic shape
    x_ = fourier_approx(maskChainCode, basicHarmonic, round(sqrt(nnz(wingMask))), 0);
    reconstructPts0 = [x_; x_(1,1) x_(1,2)]+chainBeginPt; % Make it closed contour
    reconstructPts=[reconstructPts0(:,1),size(wingreduce,1)-reconstructPts0(:,2)];

%     %Plot on the shape without ornament
%     figure,imshow(wingreduce);hold on;
%     plot(reconstShape(:,1),reconstShape(:,2),'b','LineWidth',2);
%     plot(reconstructPts(:,1), reconstructPts(:,2), 'r','LineWidth',2);
%     


%     LRside='L';
    if LRside=='L'
        lower_wing_body_join=in_key(7,:);
        keyPts=[in_key(2,:) ;in_key(6,:) ];
    elseif LRside=='R'
        lower_wing_body_join=in_key(9,:);
        keyPts=[in_key(3,:) ;in_key(5,:) ];
    end


    % dist = sqrt((reconstructPts(:,1)-lower_wing_body_join(1)).^2 + (reconstructPts(:,2)-lower_wing_body_join(2)).^2);
    % reconstructiniPt=find(dist==min(dist));
    % reconstructPts2=[reconstructPts(reconstructiniPt:length(reconstructPts),:) ; reconstructPts(1:reconstructiniPt-1,:)];
    % dist2=sqrt((reconstructPts2(:,1)-lower_wing_body_join(1)).^2 + (reconstructPts2(:,2)-lower_wing_body_join(2)).^2);

    %slope = atand((reconstructPts(:,2)-lower_wing_body_join(2)) ./ (reconstructPts(:,1)-lower_wing_body_join(1))) ;
    % figure,plot(dist);daspect([1 1 1]);
    % figure,plot(slope);daspect([1 1 1]);

    %%Determine which bounding shape is more suitable
    [sideTrix,sideTriy] = minboundtri(reconstructPts(:,1),reconstructPts(:,2));
    minBondingTriangle = polyarea(sideTrix,sideTriy);

    [rectx,recty,rectArea,rectPeri] = minboundrect(reconstructPts(:,1),reconstructPts(:,2));

    ratio2MinBondingTri=nnz(wingMask)/minBondingTriangle;
    ratio2MinBondingRect=nnz(wingMask)/rectArea;
    %[ratio2MinBondingTri ; ratio2MinBondingRect ; ratio2MinBondingRect-ratio2MinBondingTri ; (ratio2MinBondingRect-ratio2MinBondingTri)/(ratio2MinBondingRect+ratio2MinBondingTri)]
    
    if (ratio2MinBondingRect-ratio2MinBondingTri)/(ratio2MinBondingRect+ratio2MinBondingTri)>0.1 %Use rectangle module
        canRecxy=[rectx(1:end-1),recty(1:end-1)];
        if LRside=='L'
            farestPts=findClosest2Pts(canRecxy,[0,size(wingMask,1)]);
%             farestPt=findCloestPt(canRecxy,[0,size(wingMask,1)]);
%             remainRecxy=findCloestPt(canRecxy,[0,size(wingMask,1)/2]);
        elseif LRside=='R'
            farestPts=findClosest2Pts(canRecxy,flip(size(wingMask)));
%             farestPt=findCloestPt(canRecxy,flip(size(wingMask)));
%             remainRecxy=findCloestPt(canRecxy,[size(wingMask,2),size(wingMask,1)/2]);
        end
        farestPt=farestPts(farestPts(:,2)==max(farestPts(:,2)),:);
        remainRecxy=farestPts(farestPts(:,2)==min(farestPts(:,2)),:);
        
        reconstructWingTip=findCloestPt(reconstructPts,farestPt);
        reconstructSideTip=findCloestPt(reconstructPts,remainRecxy);
    else %Use triangle module
        canTrixy=[sideTrix(1:end-1),sideTriy(1:end-1)];
        [~,lowestLoc]=max(canTrixy(:,2));
        farestPt=canTrixy(lowestLoc,:);
        reconstructWingTip=findCloestPt(reconstructPts,farestPt);

        remainTrixy=canTrixy;
        remainTrixy(lowestLoc,:)=[];
        bodyJoin=findCloestPt(remainTrixy,lower_wing_body_join);
        remainTrixy(find(remainTrixy(:,1)==bodyJoin(1) & remainTrixy(:,2)==bodyJoin(2)),:)=[];
        reconstructSideTip=findCloestPt(reconstructPts,remainTrixy);
    end
    
    reconstructRegPtH=[lower_wing_body_join; reconstructWingTip; reconstructSideTip];

%     ornamentRatio=1/150;
    ornamentThreshold=nnz(wingMask)*ornamentRatio; %Define the threshold to keep for refined mask (objects will be defined as ornament among this threshold)
%     if ornamentThreshold>100, ornamentThreshold=200;, end; %Set a maximum to the threshold
    %Project the reconstruct shape on the original shape
    reconstructFFTShp = roipoly(wingMask,reconstructPts(:,1),reconstructPts(:,2));
    overlapArea = imdilate(reconstructFFTShp,strel('disk',5)) & wingMask;
    remainArea = wingMask~=overlapArea;
    keepArea = bwareafilt(remainArea,[0, ornamentThreshold]); %Keep those regions in the refine mask
    [edB,edL] = bwboundaries(bwareafilt(remainArea,[ornamentThreshold,nnz(wingMask)]),'noholes');
    
    overlapAreaRef = imdilate(reconstructFFTShp,strel('disk',10)) & wingMask;
    remainAreaRef = wingMask~=overlapAreaRef;
%     keepAreaRef = bwareafilt(remainAreaRef,[0, ornamentThreshold]);
    [refB,refL] = bwboundaries(remainAreaRef,'noholes');

    %statsed = regionprops(edL,'Area','Centroid');
    if ~isempty(refB)
        retainIdx=zeros(length(refB),3);
        for edID=1:length(refB)
            orn=refL==edID;
%             if nnz(orn)>=5
                [idx,idy]=find(orn);
%                 if ~isempty(idx)
                Lorn=bwselect(edL,idy(ceil(length(idx)/2)),idx(ceil(length(idx)/2)));
%                 figure,imshowpair(orn,Lorn);
                reducedAreaRatio=nnz(immultiply(imcomplement(orn),Lorn))/nnz(Lorn);
                if reducedAreaRatio>0.5 %If the reduced area is over 50%, the part must attached closely to the main mask but not point out
                    retainIdx(edID,:)=[mean(edL(Lorn),'all'), nnz(orn), 1];
                else
                    retainIdx(edID,:)=[mean(edL(Lorn),'all'), nnz(orn), 0];
                end
%                 else
%                     keepArea=keepArea | orn;
%                 end
%             end
        end
        retainCandidates=retainIdx(retainIdx(:,3)==1,:);
%         retainIdx2=retainIdx;
        decisionIdx=[];
        for retID=1:size(retainCandidates,1)
            subIdx=retainIdx(retainIdx(:,1)==retainCandidates(retID,1),:);
            decisionIdx0=[subIdx(1,1), subIdx(subIdx(:,2)==max(subIdx(:,2)),3)]; %Keep or not determined by the elements having the largest area
            if size(decisionIdx0,1)>1,  decisionIdx0= decisionIdx0(1,:);, end; %Added April 18, 2020 to prevent 2 candidates
            decisionIdx=[decisionIdx ; decisionIdx0];
        end
        
        if ~isempty(decisionIdx)
        retainIdxf=decisionIdx(decisionIdx(:,2)==1,1);
            if ~isempty(retainIdxf)
                retainMask=ismember(edL,retainIdxf);
                keepArea=keepArea | retainMask;
            end
        end
    end
    
    
    for edID=1:length(edB)
        obj=flip(edB{edID},2);
        wing_body_join0=findCloestPt([obj(:,1),obj(:,2)],[lower_wing_body_join(1),lower_wing_body_join(2)]);
        wing_body_seg0=findCloestPt([obj(:,1),obj(:,2)],[keyPts(2,1), keyPts(2,2)]);
        
        dist_wing_body_join=pdist([wing_body_join0;lower_wing_body_join]);
        dist_wing_body_seg=pdist([wing_body_seg0;keyPts(2,:)]);
        if dist_wing_body_join<5 | dist_wing_body_seg<5
            keepArea=keepArea | (edL==edID);
        end
    end
    
    maskArea1= overlapArea + keepArea;
    reconstructShp = roipoly(wingMask,reconstShape(:,1),reconstShape(:,2));
    remainArea2 = reconstructShp~=(reconstructShp & maskArea1);
    refineArea0= bwareafilt(imfill((maskArea1+remainArea2) & wingMask,'holes'),1);

    [Bo,~]=bwboundaries(refineArea0);
    refineAreaEdgePt=flip(Bo{1},2);
    
    wingTip=findRadiativeCloestPt(refineAreaEdgePt,reconstructWingTip);
    sideTip=findRadiativeCloestPt(refineAreaEdgePt,reconstructSideTip);
    
    %If the lower body wing corner is not on the edge, we add back the
    %body-wing join part
    [in,on] = inpolygon(lower_wing_body_join(1),lower_wing_body_join(2),refineAreaEdgePt(:,1),refineAreaEdgePt(:,2));
    if  in+on>0
        regPtH=[lower_wing_body_join; wingTip; sideTip];
        refineArea=refineArea0;
    else
        rad = min(sqrt(sum(bsxfun(@minus, [wingTip; sideTip ; keyPts(2,:)], lower_wing_body_join).^2,2)));
        wing_body_join_mask = createCirclesMask(wingMask, lower_wing_body_join, rad);
        refineArea=immultiply(refineArea0,imcomplement(wing_body_join_mask))+immultiply(wingMask,wing_body_join_mask);
        [B1,~]=bwboundaries(refineArea);
        refineAreaEdgePt2=flip(B1{1},2);
        wingTip=findRadiativeCloestPt(refineAreaEdgePt2,reconstructWingTip);
        sideTip=findRadiativeCloestPt(refineAreaEdgePt2,reconstructSideTip);
        regPtH=[lower_wing_body_join; wingTip; sideTip];
    end
% 
%     %Plot on the original shape
%     figure,imshow(wingMask);hold on;
%     plot(refineAreaEdgePt(:,1),refineAreaEdgePt(:,2),'b','LineWidth',1.5);
%     %plot(reconstructPts(:,1), reconstructPts(:,2), 'r','LineWidth',2);
%     plot(regPtH(:,1),regPtH(:,2),'gx');plot(regPtH(:,1),regPtH(:,2),'go');
% 
% 
%     figure,imshow(maskArea);hold on;
%     plot(reconstShape(:,1),reconstShape(:,2),'b','LineWidth',1);
%     plot(reconstructPts(:,1), reconstructPts(:,2), 'r','LineWidth',1);
%     plot(edgePt(:,2), edgePt(:,1), 'g','LineWidth',1);
% 
% 
%     figure,imshow(wingMask);hold on;
%     plot(sideTrix,sideTriy,'r');
%     plot(reconstructWingTip(1),reconstructWingTip(2),'rO');
%     plot(reconstructSideTip(1),reconstructSideTip(2),'rO');

end