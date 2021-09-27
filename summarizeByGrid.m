function gridSummaryMeanSE=summarizeByGrid(baseImg,wingGrids,nullArea)

gridSummaryMeanSE=cell(0,2);
if size(baseImg,3)==1
%     disp('This is a single layer gray image.');
    gridSummary=summarizeByGrid0(baseImg,wingGrids,nullArea);
    gridSummaryMeanSE{1}=gridSummary(:,:,1);
    gridSummaryMeanSE{2}=gridSummary(:,:,2);
else
%     disp('This is an RGB image.');
    gridSummaryMean=zeros(size(wingGrids,1)-1,size(wingGrids,2)-1,size(baseImg,3));
    gridSummarySE=zeros(size(wingGrids,1)-1,size(wingGrids,2)-1,size(baseImg,3));
    for layer=1:size(baseImg,3)
        imgLayer=baseImg(:,:,layer);
        gridSummaryLayer=summarizeByGrid0(imgLayer,wingGrids,nullArea);
        gridSummaryMean(:,:,layer)=gridSummaryLayer(:,:,1);
        gridSummarySE(:,:,layer)=gridSummaryLayer(:,:,2);
    end
    gridSummaryMeanSE{1}=gridSummaryMean;
    gridSummaryMeanSE{2}=gridSummarySE;
end



%%
    function gridSummaryMeanSE0=summarizeByGrid0(baseImg,wingGrids,nullArea)
        gridSummaryMeanSE0=zeros(size(wingGrids,1)-1,size(wingGrids,2)-1,2);
        for i=1:size(wingGrids,1)-1
            for j=1:size(wingGrids,2)-1
                %disp([num2str(i),' - ',num2str(j)]);
                cornerData=reshape([wingGrids(i,j,:) ; wingGrids(i,j+1,:) ; wingGrids(i+1,j,:) ; wingGrids(i+1,j+1,:)],[],2);
                [x2,y2]=orderPtsPoly(cornerData(:,1), cornerData(:,2));
                cropGrid = roipoly(baseImg,x2,y2);
                if nnz(cropGrid)>0
                    gridData0=baseImg(cropGrid==1);
                    nullData=nullArea(cropGrid==1);
                    gridData=gridData0(nullData==0); %Remove those irrelavant data
                    gridNullOrNot=mode(nullData);
                    if gridNullOrNot==0
                        avgvalue=mean(gridData);
                        se=std(gridData) / sqrt( length(gridData));
                    else
                        avgvalue=-9999;
                        se=-9999;
                    end
                else
                    avgvalue=-9999;
                    se=-9999;
                end
                gridSummaryMeanSE0(i,j,1)=avgvalue;
                gridSummaryMeanSE0(i,j,2)=se;
            end
        end
    end
end