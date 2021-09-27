function antImgf=antennaImg2(antPlotMat_dorsal, bandGatherList, antImgSize, avgBarWidth)
%     antImgSize=[1000,300];
%     bandGatherList={3,[6,7,8],1,2,[15,16,17],[18,19,20]};
% avgBarWidth=10;
if ~isempty(antPlotMat_dorsal) && nanmean(antPlotMat_dorsal, 'all')~=0 && nanmean(~isnan(antPlotMat_dorsal), 'all')~=0
    antImg=[];
    midBars=[];
    minBars=[];
    maxBars=[];
    for colID=1:length(bandGatherList)
        colIdx0=bandGatherList{colID};
        if length(colIdx0)<3
            colIdx=zeros(1,3)+colIdx0;
        else
            colIdx=colIdx0;
        end
        
        if length(colIdx0)<3
        bandRef1=[];
            for colid=1:3
                bandRef0=imresize(antPlotMat_dorsal(:,colIdx(colid)),[antImgSize(1),round(antImgSize(2)/6)]);
                bandRef1=cat(3,bandRef1,imadjust(bandRef0));
            end
            bandRef=bandRef1;
            midBar=zeros(avgBarWidth,size(bandRef,2),3)+nanmedian(antPlotMat_dorsal(:,colIdx(1)));
            minBar=zeros(avgBarWidth,size(bandRef,2),3)+nanmin(antPlotMat_dorsal(:,colIdx(1)));
            maxBar=zeros(avgBarWidth,size(bandRef,2),3)+nanmax(antPlotMat_dorsal(:,colIdx(1)));
        else
            bandRef1=[];
            minmaxList=zeros(3,2);
            for colid=1:3
                bandRef0=imresize(antPlotMat_dorsal(:,colIdx(colid)),[antImgSize(1),round(antImgSize(2)/6)]);
                bandRef1=cat(3,bandRef1,bandRef0);
                minmaxList(colid,:)=[min(bandRef0,[],'all'), max(bandRef0,[],'all')];
            end
            bandRef=imadjust(bandRef1,[0 0 0 ; max(minmaxList(:,2)), max(minmaxList(:,2)), max(minmaxList(:,2))],[]);
            midBarR=zeros(avgBarWidth,size(bandRef,2))+nanmedian(bandRef1(:,:,1), 'all');
            midBarG=zeros(avgBarWidth,size(bandRef,2))+nanmedian(bandRef1(:,:,2), 'all');
            midBarB=zeros(avgBarWidth,size(bandRef,2))+nanmedian(bandRef1(:,:,3), 'all');
            midBar=cat(3, midBarR, midBarG, midBarB);
            
            minBarR=zeros(avgBarWidth,size(bandRef,2))+nanmin(bandRef1(:,:,1));
            minBarG=zeros(avgBarWidth,size(bandRef,2))+nanmin(bandRef1(:,:,2));
            minBarB=zeros(avgBarWidth,size(bandRef,2))+nanmin(bandRef1(:,:,3));
            minBar=cat(3, minBarR, minBarG, minBarB);
            
            maxBarR=zeros(avgBarWidth,size(bandRef,2))+nanmax(bandRef1(:,:,1));
            maxBarG=zeros(avgBarWidth,size(bandRef,2))+nanmax(bandRef1(:,:,2));
            maxBarB=zeros(avgBarWidth,size(bandRef,2))+nanmax(bandRef1(:,:,3));
            maxBar=cat(3, maxBarR, maxBarG, maxBarB);
        end
        antImg=cat(2, antImg, bandRef);
        midBars=cat(2, midBars, midBar);
        minBars=cat(2, minBars, minBar);
        maxBars=cat(2, maxBars, maxBar);
%         disp(num2str(size(imgR,2)));
    end
else
    antImg=zeros(antImgSize(1),antImgSize(2),3);
    minBars=zeros(avgBarWidth,antImgSize(2),3);
    midBars=zeros(avgBarWidth,antImgSize(2),3);
    maxBars=zeros(avgBarWidth,antImgSize(2),3);
end

 antImgf=cat(1, antImg, maxBars, midBars, minBars);
end