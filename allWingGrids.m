function gridsParameter=allWingGrids(wingLF,wingRF,wingLH,wingRH,in_key,numberOfIntervalDegree)
%ventral image
% wingLF=potentialWingMask_LF_ventral;
% wingRF=potentialWingMask_RF_ventral;
% wingLH=ventral_seg_flip==3;
% wingRH=ventral_seg_flip==4;

ornamentRatio=1/350; %the ratio of object's area to mask area exceeding this value will be defined as an ornament
[refineAreaLH,regPtLH,reconstructRegPtLH]=refineHindWing2(wingLH,in_key,'L', ornamentRatio); %Left and right indicated by how the image currently looks
[refineAreaRH,regPtRH,reconstructRegPtRH]=refineHindWing2(wingRH,in_key,'R', ornamentRatio); %Left and right indicated by how the image currently looks

% figure,imshowpair(wingLH,refineAreaLH);hold on;
% plot(reconstructRegPtLH(:,1),reconstructRegPtLH(:,2),'rO');
% figure,imshowpair(wingRH,refineAreaRH);hold on;
% plot(reconstructRegPtRH(:,1),reconstructRegPtRH(:,2),'rO');

%generate grids for forewings
[seg4PtsLF,wingGridsLF ]=foreWingGrids(wingLF,in_key,'L',numberOfIntervalDegree);
[seg4PtsRF,wingGridsRF ]=foreWingGrids(wingRF,in_key,'R',numberOfIntervalDegree);

%generate grids for hindwings
[outSeg4PtsLH ,wingGridsLH ]=hindWingGrids(refineAreaLH,in_key,regPtLH,'L',numberOfIntervalDegree);
[outSeg4PtsRH ,wingGridsRH ]=hindWingGrids(refineAreaRH,in_key,regPtRH,'R',numberOfIntervalDegree);


gridsParameter=cell(0,6);
gridsParameter{1}={seg4PtsLF,wingGridsLF };
gridsParameter{2}={seg4PtsRF,wingGridsRF};
gridsParameter{3}={outSeg4PtsLH ,wingGridsLH};
gridsParameter{4}={outSeg4PtsRH ,wingGridsRH};
gridsParameter{5}={refineAreaLH,regPtLH,reconstructRegPtLH};
gridsParameter{6}={refineAreaRH,regPtRH,reconstructRegPtRH};
gridsParameter{7}=wingLH;
gridsParameter{8}=wingRH;


% %Plot refined regions
% figure,imshowpair(wingLH+wingRH,refineAreaLH+refineAreaRH);hold on;
% plot(regPtLH(:,1),regPtLH(:,2),'rX');plot(regPtLH(:,1),regPtLH(:,2),'yO');
% plot(regPtRH(:,1),regPtRH(:,2),'rX');plot(regPtRH(:,1),regPtRH(:,2),'yO');
% 
% %plot all grids on forewings
% figure,imshow(wingLF+wingRF);hold on;
% for i=2:size(wingGridsLF,1)-1
%     gridPlot=reshape(wingGridsLF(i,:,:),[],2);
%     plot(gridPlot(:,1),gridPlot(:,2),'r');
% end
% for j=2:size(wingGridsLF,2)-1
%     gridPlot=reshape(wingGridsLF(:,j,:),[],2);
%     plot(gridPlot(:,1),gridPlot(:,2),'r');
% end
% for i=2:size(wingGridsRF,1)-1
%     gridPlot=reshape(wingGridsRF(i,:,:),[],2);
%     plot(gridPlot(:,1),gridPlot(:,2),'r');
% end
% for j=2:size(wingGridsRF,2)-1
%     gridPlot=reshape(wingGridsRF(:,j,:),[],2);
%     plot(gridPlot(:,1),gridPlot(:,2),'r');
% end
% 
% %plot all grids on hindwings
% figure,imshow(wingLH+wingRH);hold on;
% for i=2:size(wingGridsLH,1)-1
%     gridPlot=reshape(wingGridsLH(i,:,:),[],2);
%     plot(gridPlot(:,1),gridPlot(:,2),'r');
% end
% for j=2:size(wingGridsLH,2)-1
%     gridPlot=reshape(wingGridsLH(:,j,:),[],2);
%     plot(gridPlot(:,1),gridPlot(:,2),'r');
% end
% for i=2:size(wingGridsRH,1)-1
%     gridPlot=reshape(wingGridsRH(i,:,:),[],2);
%     plot(gridPlot(:,1),gridPlot(:,2),'r');
% end
% for j=2:size(wingGridsRH,2)-1
%     gridPlot=reshape(wingGridsRH(:,j,:),[],2);
%     plot(gridPlot(:,1),gridPlot(:,2),'r');
% end


end