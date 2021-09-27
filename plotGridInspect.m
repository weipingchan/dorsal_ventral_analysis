function plotGridInspect(wingMaskLF,wingMaskRF,wingMaskLH,wingMaskRH,nullArea_L,nullArea_R,gridsParameter)
% wingMaskLF=potentialWingMask_LF_ventral;
% wingMaskRF=potentialWingMask_RF_ventral;
% wingMaskLH=ventral_seg==3;
% wingMaskRH=ventral_seg==4;
% nullArea_L=nullArea_LF_ventral;
% nullArea_R=nullArea_RF_ventral;
% gridsParameter=gridsParameter_ventral;

plotbaseImg=wingMaskLF+wingMaskRF+wingMaskLH+wingMaskRH;
plotrefineImg=nullArea_L+nullArea_R+gridsParameter{5}{1}+gridsParameter{6}{1};
seg4PtsAll=[gridsParameter{1}{1} ; gridsParameter{2}{1} ; gridsParameter{3}{1} ; gridsParameter{4}{1}];

wingGridsLF=gridsParameter{1}{2};
wingGridsRF=gridsParameter{2}{2};
wingGridsLH=gridsParameter{3}{2};
wingGridsRH=gridsParameter{4}{2};

imshowpair(plotbaseImg,plotrefineImg);hold on;
plot(seg4PtsAll(:,1),seg4PtsAll(:,2),'r+'); plot(seg4PtsAll(:,1),seg4PtsAll(:,2),'ro');
%plot all grids on Fore wings
for i=2:size(wingGridsLF,1)-1
    gridPlot=reshape(wingGridsLF(i,:,:),[],2);
    plot(gridPlot(:,1),gridPlot(:,2),'r');
end
for j=2:size(wingGridsLF,2)-1
    gridPlot=reshape(wingGridsLF(:,j,:),[],2);
    plot(gridPlot(:,1),gridPlot(:,2),'r');
end
for i=2:size(wingGridsRF,1)-1
    gridPlot=reshape(wingGridsRF(i,:,:),[],2);
    plot(gridPlot(:,1),gridPlot(:,2),'r');
end
for j=2:size(wingGridsRF,2)-1
    gridPlot=reshape(wingGridsRF(:,j,:),[],2);
    plot(gridPlot(:,1),gridPlot(:,2),'r');
end
%plot all grids on Hind wings
for i=2:size(wingGridsLH,1)-1
    gridPlot=reshape(wingGridsLH(i,:,:),[],2);
    plot(gridPlot(:,1),gridPlot(:,2),'b');
end
for j=2:size(wingGridsLH,2)-1
    gridPlot=reshape(wingGridsLH(:,j,:),[],2);
    plot(gridPlot(:,1),gridPlot(:,2),'b');
end
for i=2:size(wingGridsRH,1)-1
    gridPlot=reshape(wingGridsRH(i,:,:),[],2);
    plot(gridPlot(:,1),gridPlot(:,2),'b');
end
for j=2:size(wingGridsRH,2)-1
    gridPlot=reshape(wingGridsRH(:,j,:),[],2);
    plot(gridPlot(:,1),gridPlot(:,2),'b');
end
hold off;
end