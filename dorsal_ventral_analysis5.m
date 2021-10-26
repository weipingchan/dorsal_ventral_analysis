function dorsal_ventral_analysis5(morph_mat_directory,spp_mat_directory,manual_grids_directory,Code_directory,Result_directory,barcodein,numberOfIntervalDegree)
% Turn off this warning "Warning: Image is too big to fit on screen; displaying at 33% "
% To set the warning state, you must first know the message identifier for the one warning you want to enable. 
warning('off', 'Images:initSize:adjustingMag');
addpath(genpath(Code_directory)) %Add the library to the path

vdlist={'dorsal','ventral'};

%Convert double quotes to single quotes for the Matlab version prior to
%2017
if size(morph_mat_directory,2)==1 morph_mat_directory=morph_mat_directory{1};, end;
if size(spp_mat_directory,2)==1 spp_mat_directory=spp_mat_directory{1};, end;
if size(Code_directory,2)==1 Code_directory=Code_directory{1};, end;
if size(Result_directory,2)==1 Result_directory=Result_directory{1};, end;
if size(barcodein,2)==1 barcodein=barcodein{1};, end;
% if size(numberOfIntervalDegree,2)==1 numberOfIntervalDegree=numberOfIntervalDegree{1};, end;
if ~isnumeric(numberOfIntervalDegree) numberOfIntervalDegree=str2num(numberOfIntervalDegree);, end;


disp('Start to create / find corresponding folders.');
%Create result directory
if ~exist(fullfile(Result_directory,'dorsal_ventral_map'), 'dir')
    mkdir(fullfile(Result_directory,'dorsal_ventral_map'));
end

subFolderList={'inspect_imgs','wing_matrix_visualization','spp_wing_parameters'};
subsubFolderList={'740','940','UV','UVF','F','white','whitePo1','whitePo2','FinRGB','PolDiff','antenna'};

for fold=1:length(subFolderList)
    if ~exist(fullfile(Result_directory,'dorsal_ventral_map',subFolderList{fold}), 'dir')
        mkdir(fullfile(Result_directory,'dorsal_ventral_map',subFolderList{fold}));
    end
    disp(['corresponding folder ', subFolderList{fold}, ' is created / found.']);
end

for fold=1:length(subsubFolderList)
    if ~exist(fullfile(Result_directory,'dorsal_ventral_map',subFolderList{2},['b_',subsubFolderList{fold}]), 'dir')
        mkdir(fullfile(Result_directory,'dorsal_ventral_map',subFolderList{2},['b_',subsubFolderList{fold}]));
    end
    disp(['corresponding folder ', subsubFolderList{fold}, ' is created / found.']);
end


%%
all_morph_data = dir(fullfile(morph_mat_directory,'*morph-seg.mat')); %read only one layer of directory; '**/*morph-seg.mat': read all subfolders
all_morph_name0=struct2dataset(all_morph_data);
all_morph_name=all_morph_name0(:,1).name;
fileID=find(contains(all_morph_name,[barcodein,'_']));
both_sides_morph=cell(0,2);
if length(fileID)==2
    for matinID=1:length(fileID)
        matindir=all_morph_data(fileID(matinID)).folder;
        matinname=all_morph_data(fileID(matinID)).name;
        pathdirs=strsplit(matindir,filesep);
        lowestdir=pathdirs{end};
        [barcode, side, flag]=file_name_decoder(matinname);
        matin=fullfile(matindir,matinname);
        sppmat0=load(matin);
        fieldName=cell2mat(fieldnames(sppmat0));
        sppmat=sppmat0.(fieldName);
        clear sppmat0;
        disp(['[',matinname,'] in {',lowestdir,'} has been read into memory']);

        both_sides_morph{side}={sppmat,barcode,side,flag};
    end
else
    disp('cannot find both sides of morpho-seg mat files');
    matinname=all_morph_data(fileID).name;
    disp(['[',matinname,']']);
end

%%
%derive the matchable side of images
dorsal_seg=both_sides_morph{1}{1}{13};
ventral_seg=both_sides_morph{2}{1}{13};
ventral_seg_flip = flip(ventral_seg ,2); 

keyPts=cell(0,2);
for sideID=1:2
    refPts=both_sides_morph{sideID}{1}{6};
    segPts=refPts(1:6,:);
    realCen=both_sides_morph{sideID}{1}{2};
    tipPts=both_sides_morph{sideID}{1}{5};
    keyPts{sideID}=[refPts;tipPts];
end

dorsal_key=keyPts{1};
ventral_key=keyPts{2};
ventral_key_flip0=ventral_key;
ventral_key_flip0(:,1)=size(ventral_seg,2)+1-ventral_key(:,1);
ventral_key_flip=[ventral_key_flip0(4,:);ventral_key_flip0(3,:);ventral_key_flip0(2,:);ventral_key_flip0(1,:);ventral_key_flip0(6,:);ventral_key_flip0(5,:);ventral_key_flip0(9,:);ventral_key_flip0(8,:);ventral_key_flip0(7,:);ventral_key_flip0(11,:);ventral_key_flip0(10,:)];

dorsal_scale=both_sides_morph{1}{1}{12};
ventral_scale=both_sides_morph{2}{1}{12};
%Create potential intact wing region based on dorsal-ventral mapping
disp('Begin to create potential wing area based on dosal-ventral sides mapping');

%Check all wings
%If wing is missing, provide a fake one
dorsal_wingLF=checkWingMask(dorsal_seg==1, dorsal_key, 'LF');
dorsal_wingRF=checkWingMask(dorsal_seg==2, dorsal_key, 'RF');
dorsal_wingLH=checkWingMask(dorsal_seg==3, dorsal_key, 'LH');
dorsal_wingRH=checkWingMask(dorsal_seg==4, dorsal_key, 'RH');

ventral_flip_wingLF=checkWingMask(ventral_seg_flip==1, ventral_key_flip, 'LF');
ventral_flip_wingRF=checkWingMask(ventral_seg_flip==2, ventral_key_flip, 'RF');
ventral_flip_wingLH=checkWingMask(ventral_seg_flip==3, ventral_key_flip, 'LH');
ventral_flip_wingRH=checkWingMask(ventral_seg_flip==4, ventral_key_flip, 'RH');

%Use dorsal image
segLineLF=extractSegLineF(dorsal_wingLF,dorsal_key,'L');
segLineRF=extractSegLineF(dorsal_wingRF,dorsal_key,'R');

%Use ventral-flip image
segLineLH=extractSegLineH(ventral_flip_wingLH,ventral_key_flip,'L');
segLineRH=extractSegLineH(ventral_flip_wingRH,ventral_key_flip,'R');

disp('The fore-hindwing joint lines are extracted from dosal and ventral sides');
% figure,imshow(wingLF);hold on;
% plot(segLineLF(:,1),segLineLF(:,2),'r.');
% 
% figure,imshow(wingRF);hold on;
% plot(segLineRF(:,1),segLineRF(:,2),'r.');
% 
% figure,imshow(wingLH);hold on;
% plot(segLineLH(:,1),segLineLH(:,2),'r.');
% 
% figure,imshow(wingRH);hold on;
% plot(segLineRH(:,1),segLineRH(:,2),'r.');

%Map ventral-side hindwing on dorsal side hindwing
segLineLH_dorsal=projectSegLine(dorsal_wingLH,dorsal_key,segLineLH,'L');
segLineRH_dorsal=projectSegLine(dorsal_wingRH,dorsal_key,segLineRH,'R');
segLineLF_ventral=projectSegLine(ventral_flip_wingLF,ventral_key_flip,segLineLF,'L');
segLineRF_ventral=projectSegLine(ventral_flip_wingRF,ventral_key_flip,segLineRF,'R');

% figure,imshow((ventral_seg_flip==1)+(ventral_seg_flip==2));hold on;
% plot(segLineLF_ventral(:,1),segLineLF_ventral(:,2),'r.');
% plot(segLineRF_ventral(:,1),segLineRF_ventral(:,2),'r.');
% 
% figure,imshow((dorsal_seg==3)+(dorsal_seg==4));hold on;
% plot(segLineLH_dorsal(:,1),segLineLH_dorsal(:,2),'r.');
% plot(segLineRH_dorsal(:,1),segLineRH_dorsal(:,2),'r.');

%Generate potential mask
[potentialWingMask_LH_dorsal,nullArea_LH_dorsal]=getPotentialWingMask2(dorsal_wingLH,segLineLH_dorsal,'LH');
[potentialWingMask_RH_dorsal,nullArea_RH_dorsal]=getPotentialWingMask2(dorsal_wingRH,segLineRH_dorsal,'RH');
[potentialWingMask_LF_ventral,nullArea_LF_ventral]=getPotentialWingMask2(ventral_flip_wingLF,segLineLF_ventral,'LF');
[potentialWingMask_RF_ventral,nullArea_RF_ventral]=getPotentialWingMask2(ventral_flip_wingRF,segLineRF_ventral,'RF');

% figure,imshowpair(potentialWingMask_LH_dorsal+potentialWingMask_RH_dorsal,nullArea_LH_dorsal+nullArea_RH_dorsal);
% figure,imshowpair(potentialWingMask_LF_ventral+potentialWingMask_RF_ventral,nullArea_LF_ventral+nullArea_RF_ventral);
disp('The potential wing area has been estimated');

%%
% numberOfIntervalDegree=5; %minimum is 3
disp(['The resolution of wing matrix is: ' num2str(2^numberOfIntervalDegree),' x ',num2str(2^numberOfIntervalDegree)]);

%Search for pre-defined grid files
manualGridin=[];
manual_grids_data = dir(fullfile(manual_grids_directory,'*d-v_manual_grids.mat'));
for gridinID=1:length(manual_grids_data)
    gridinname=manual_grids_data(gridinID).name;
    gridBarcode0=strsplit(gridinname,'_');
    gridBarcode=gridBarcode0{1};
    res0=strsplit(gridBarcode0{2},'x');
    res=str2num(res0{end});
    
    if strcmp(gridBarcode,barcodein) && res==2^numberOfIntervalDegree
        manualGridin=gridinname;
        break
    end
end

if ~isempty(manualGridin)
        disp('Find corresponding grid files. Begin to use pre-defined grids.');
        gridmat0=load(fullfile(manual_grids_directory, manualGridin));
        fieldName0=cell2mat(fieldnames(gridmat0));
        gridmat=gridmat0.(fieldName0);
        disp(['[',manualGridin,'] has been read into memory']);
        gridsParameter_dorsal=gridmat{1};
        disp('Grids of dorsal wing are derived');
        gridsParameter_ventral=gridmat{2};
        disp('Grids of ventral wing are derived');
else
    disp('Did not find corresponding grid files. Begin to generate grids.');
    %dorsal image
    gridsParameter_dorsal=allWingGrids(dorsal_wingLF,dorsal_wingRF,potentialWingMask_LH_dorsal,potentialWingMask_RH_dorsal,dorsal_key,numberOfIntervalDegree);
    disp('Grids of dorsal wing have been generated');

    %ventral image
    gridsParameter_ventral=allWingGrids(potentialWingMask_LF_ventral,potentialWingMask_RF_ventral,ventral_flip_wingLH,ventral_flip_wingRH,ventral_key_flip,numberOfIntervalDegree);
    disp('Grids of ventral wing have been generated');
end
% {1}={seg4PtsLF,wingGridsLF }; %4 key points & grids
% {2}={seg4PtsRF,wingGridsRF}; %4 key points & grids
% {3}={outSeg4PtsLH ,wingGridsLH}; %4 key points & grids
% {4}={outSeg4PtsRH ,wingGridsRH}; %4 key points & grids
% {5}={refineAreaLH,regPtLH,reconstructRegPtLH}; %mask without ornament, 3
% key points & key points of the mask without ornament
% {6}={refineAreaRH,regPtRH,reconstructRegPtRH}; %mask without ornament, 3
% key points & key points of the mask without ornament
% {7}=original mask of LH
% {8}=original mask of RH

%Save image for inspection
%Save dorsal side
inspvisoutname=fullfile(Result_directory,'dorsal_ventral_map',subFolderList{1},[barcode,'_',vdlist{1},flag,'_keys_grids_res-',num2str(2^numberOfIntervalDegree),'x',num2str(2^numberOfIntervalDegree),'.jpg']);
figinsp=figure('visible', 'off');
plotGridInspect(dorsal_wingLF,dorsal_wingRF,potentialWingMask_LH_dorsal,potentialWingMask_RH_dorsal,nullArea_LH_dorsal,nullArea_RH_dorsal,gridsParameter_dorsal);
%saveas(figmask, maskoutname);
export_fig(figinsp,inspvisoutname, '-jpg','-r200');
close(figinsp);

%Save ventral side
inspvisoutname=fullfile(Result_directory,'dorsal_ventral_map',subFolderList{1},[barcode,'_',vdlist{2},flag,'_keys_grids_res-',num2str(2^numberOfIntervalDegree),'x',num2str(2^numberOfIntervalDegree),'.jpg']);
figinsp=figure('visible', 'off');
plotGridInspect(potentialWingMask_LF_ventral,potentialWingMask_RF_ventral,ventral_flip_wingLH,ventral_flip_wingRH,nullArea_LF_ventral,nullArea_RF_ventral,gridsParameter_ventral);
%saveas(figmask, maskoutname);
export_fig(figinsp,inspvisoutname, '-jpg','-r200');
close(figinsp);
disp('Two images showing keys and grids of both sides of a specimen have been saved.');

%%
%Tail analysis
disp('Start to analyze tails.');
tailMinArea=20; %The value here determines how small of a region can be defined as a tail
 tailinfo_dorsal_L=tail_module2(gridsParameter_dorsal, dorsal_scale, tailMinArea, 'L');
 tailinfo_dorsal_R=tail_module2(gridsParameter_dorsal, dorsal_scale, tailMinArea, 'R');
 tailinfo_ventral_L=tail_module2(gridsParameter_ventral, ventral_scale, tailMinArea, 'L');
 tailinfo_ventral_R=tail_module2(gridsParameter_ventral, ventral_scale, tailMinArea, 'R');
 
 %(-9999 represents no tail)
 
 tail_morph={tailinfo_dorsal_L, tailinfo_ventral_L, tailinfo_dorsal_R, tailinfo_ventral_R};
disp('Tails done.');
%%
disp('Start to analyze antennae and body.');
%Antennae data, Body data
%The body length and width (in cm)
%Antennae length, width, bulb width, degree of curvature (all in cm); first row is left one, second is right one.
body_morph=[both_sides_morph{1}{1}{14};both_sides_morph{2}{1}{14}];
dorsal_ant=both_sides_morph{1}{1}{15};
ventral_ant=both_sides_morph{2}{1}{15};
body_ant_morph={body_morph,[dorsal_ant(1,:);ventral_ant(1,:)],[dorsal_ant(2,:);ventral_ant(2,:)]};
%{1}=body length and width (in cm) (dorsal and ventral side in two rows)
%{2}= Antenna length, width, bulb width, degree of curvature (all in mm); Left antenna; first row is dorsal side, second is ventral side.
%{3}= Antenna length, width, bulb width, degree of curvature (all in mm); Right antenna; first row is dorsal side, second is ventral side.

%Antennae mask
%Use dorsal image
ants_d=dorsal_seg==6;
%Use ventral-flip image
ants_v=ventral_seg_flip==6;
disp('Antennae done.');
%%
try
    try
        %%
        %Grid reflectance part
        sppBothSidesGridSummary=cell(0,2);
        matInNames=cell(0,2);
        %Read images and summarize
        vdlist={'dorsal','ventral'};
        for sideID=1:2
            %Read the corresponding multi-spectral matrix
            %spp_data_list = dir(fullfile(spp_mat_directory,[barcode,'_',vdlist{sideID},'*AllBandsMask.mat']));
            spp_data_list = dir(fullfile(spp_mat_directory,[both_sides_morph{sideID}{2},'_',vdlist{both_sides_morph{sideID}{3}},both_sides_morph{sideID}{4},'_AllBandsMask.mat']));
            if isempty(spp_data_list) %If the specific file cannot be found, use a similar one with different flag
                spp_data_list = dir(fullfile(spp_mat_directory,[both_sides_morph{sideID}{2},'_',vdlist{both_sides_morph{sideID}{3}},'*_AllBandsMask.mat']));
            end
            %generate a unique barcode list
            matinname=spp_data_list.name;
            [barcode, side, flag]=file_name_decoder(matinname);
            matin=fullfile(spp_mat_directory,matinname);
            sppmat0=load(matin);
            fieldName=cell2mat(fieldnames(sppmat0));
            sppmat=sppmat0.(fieldName);
            clear sppmat0;
            matInNames{sideID}=matinname;
            disp([matinname,' has been read into system and is ready for processing.']);
            
            disp('Begin to summarize multi-spectral reflectance into grids');
            sppAllGridSummary=cell(0,10);
            for bandID=1:10
                if side==1
                    baseImg=sppmat{bandID};
                elseif side==2
                    baseImg=flip(sppmat{bandID},2);    
                end

                if side==1
                    wingGridsLF=gridsParameter_dorsal{1}{2};
                    wingGridsRF=gridsParameter_dorsal{2}{2};
                    wingGridsLH=gridsParameter_dorsal{3}{2};
                    wingGridsRH=gridsParameter_dorsal{4}{2};

                    nullAreaLF=zeros(size(baseImg,1),size(baseImg,2));
                    nullAreaRF=zeros(size(baseImg,1),size(baseImg,2));
                    nullAreaLH=nullArea_LH_dorsal;
                    nullAreaRH=nullArea_RH_dorsal;
                elseif side==2
                    wingGridsLF=gridsParameter_ventral{1}{2};
                    wingGridsRF=gridsParameter_ventral{2}{2};
                    wingGridsLH=gridsParameter_ventral{3}{2};
                    wingGridsRH=gridsParameter_ventral{4}{2};

                    nullAreaLF=nullArea_LF_ventral;
                    nullAreaRF=nullArea_RF_ventral;
                    nullAreaLH=zeros(size(baseImg,1),size(baseImg,2));
                    nullAreaRH=zeros(size(baseImg,1),size(baseImg,2));
                end

                %Grid summary for forewings
                gridSummaryMeanSELF=summarizeByGrid(baseImg,wingGridsLF,nullAreaLF);
                gridSummaryMeanSERF=summarizeByGrid(baseImg,wingGridsRF,nullAreaRF);
                gridSummaryMeanSELH=summarizeByGrid(baseImg,wingGridsLH,nullAreaLH);
                gridSummaryMeanSERH=summarizeByGrid(baseImg,wingGridsRH,nullAreaRH);

                sppAllGridSummary{bandID}={gridSummaryMeanSELF,gridSummaryMeanSERF,gridSummaryMeanSELH,gridSummaryMeanSERH};

                %Save visualizations
                if numberOfIntervalDegree<=5
                    gapPixel=2;
                else
                    gapPixel=5;
                end
                comMatrixMean=deriveWingGridMat(gapPixel,gridSummaryMeanSELF{1},gridSummaryMeanSERF{1},gridSummaryMeanSELH{1},gridSummaryMeanSERH{1});
                comMatrixSE=deriveWingGridMat(gapPixel,gridSummaryMeanSELF{2},gridSummaryMeanSERF{2},gridSummaryMeanSELH{2},gridSummaryMeanSERH{2});

                %Save mean
                inspvisoutname=fullfile(Result_directory,'dorsal_ventral_map',subFolderList{2},['b_',subsubFolderList{bandID}],[barcode,'_',vdlist{sideID},flag,'_',subsubFolderList{bandID},'_res-',num2str(2^numberOfIntervalDegree),'x',num2str(2^numberOfIntervalDegree),'_GridSum-mean.jpg']);
                figinsp=figure('visible', 'off');
                if size(comMatrixMean,3)==1
                    imshow(imadjust(comMatrixMean));
                else
                    imshow(cat(3,imadjust(comMatrixMean(:,:,1)),imadjust(comMatrixMean(:,:,2)),imadjust(comMatrixMean(:,:,3))));
                end
                %saveas(figmask, maskoutname);
                export_fig(figinsp,inspvisoutname, '-jpg','-r200');
                close(figinsp);

                %Save SE
                inspvisoutname=fullfile(Result_directory,'dorsal_ventral_map',subFolderList{2},['b_',subsubFolderList{bandID}],[barcode,'_',vdlist{sideID},flag,'_',subsubFolderList{bandID},'_res-',num2str(2^numberOfIntervalDegree),'x',num2str(2^numberOfIntervalDegree),'_GridSum-se.jpg']);
                figinsp=figure('visible', 'off');
                if size(comMatrixSE,3)==1
                    imshow(imadjust(comMatrixSE));
                else
                    imshow(cat(3,imadjust(comMatrixSE(:,:,1)),imadjust(comMatrixSE(:,:,2)),imadjust(comMatrixSE(:,:,3))));
                end
                %saveas(figmask, maskoutname);
                export_fig(figinsp,inspvisoutname, '-jpg','-r200');
                close(figinsp);

                disp(['Band ',num2str(bandID),' out of ',num2str(10),' is done.']);
            end
            sppBothSidesGridSummary{sideID}=sppAllGridSummary;
            disp(['################################']);
            disp([vdlist{sideID},' side of ',barcode,' is done.' ]);
            disp(['################################']);
        end
        procedureFlag=0;
    catch
        disp('Cannot process wing reflectance summary');
        procedureFlag=1;
    end
    
    try
        %%
        %Antennae extraction
        sppBothSidesAntSummary=extractAntennaReflectance(spp_mat_directory,both_sides_morph, ants_d, ants_v, dorsal_key, ventral_key_flip, dorsal_scale, ventral_scale);

        %Antennae normalization
        disp('Start to do antennae reflectance summary.');
        outPutLevel=100;
        [antennaLRref_rescale, antennaLRrefmm]=antennaNormalization2(sppBothSidesAntSummary, outPutLevel);
        disp('Antennae reflectance summary done.');
        disp('Begin to plot antennae reflectance summary.');
        %Plot antennae
        antPlotData_L=antennaLRref_rescale{1};
        antPlotData_R=antennaLRref_rescale{2};
        try
            antPlotMat_dorsal=nanmean(cat(3,antPlotData_L(:,2:end,1), antPlotData_R(:,2:end,1)),3);
            antFlag_d=0;
            if isnan(mean(antPlotMat_dorsal,'all'))
                antPlotMat_dorsal=[];
                antFlag_d=1;
            end
        catch
            antPlotMat_dorsal=[];
            antFlag_d=1;
        end
        try
            antPlotMat_ventral=nanmean(cat(3,antPlotData_L(:,2:end,2), antPlotData_R(:,2:end,2)),3);
             antFlag_v=0;
             if isnan(mean(antPlotMat_ventral,'all'))
                antPlotMat_ventral=[];
                antFlag_v=1;
            end
        catch
            antPlotMat_ventral=[];
            antFlag_v=1;
        end

        antImgSize=[500,300];
        avgBarWidth=5;
        bandGatherList={3,[6,7,8],1,2,[15,16,17],[18,19,20]};
    %     antImg_d=antennaImg(antPlotMat_dorsal, bandGatherList, antImgSize);
    %     antImg_v=antennaImg(antPlotMat_ventral, bandGatherList, antImgSize);
        antImg_d=antennaImg2(antPlotMat_dorsal, bandGatherList, antImgSize, avgBarWidth); %Provide also original mean, min, and max reflectance
        antImg_v=antennaImg2(antPlotMat_ventral, bandGatherList, antImgSize, avgBarWidth); %Provide also original mean, min, and max reflectance

        antImg=cat(1,flip(antImg_d), antImg_v); %Dorsal side at top, ventral side at bottom, tip in the middle

        %Save antennae mean
        antvisoutname=fullfile(Result_directory,'dorsal_ventral_map',subFolderList{2},['b_',subsubFolderList{11}],[barcode,'_',subsubFolderList{11},'_res-',num2str(outPutLevel),'_dv_allBands-mean_rescale.jpg']);
        figinsp=figure('visible', 'off');
        imshow(antImg);hold on;
        plot([0, antImgSize(2)],[0+3*avgBarWidth, 0+3*avgBarWidth],'black', 'LineWidth',2);
        plot([0, antImgSize(2)],[2*antImgSize(1)+3*avgBarWidth, 2*antImgSize(1)+3*avgBarWidth],'black', 'LineWidth',2);
        if antFlag_d==1
            plot([0, antImgSize(2)],[0,antImgSize(1) ],'white');
        end
        if antFlag_v==1
            plot([0, antImgSize(2)],[antImgSize(1), 2*antImgSize(1) ],'white');
        end
        export_fig(figinsp,antvisoutname, '-jpg','-r200');
        close(figinsp);
        disp('Antennae reflectance summary has been saved');
    catch
        disp('Cannot process antennae reflectance summary');
        procedureFlag=procedureFlag+2;
    end
    %%
    processingFlagList={'', '-no_wingReflect', '-no_antennaReflect', '-no_reflectance'};
    
    disp('Begin to save all summary matrices.');
    allResult=cell(0,7);
    allResult{1}={gridsParameter_dorsal,gridsParameter_ventral}; %All parameters and grids
    if procedureFlag~=1 && procedureFlag~=3
        allResult{2}=sppBothSidesGridSummary; %All reflectance grids summary
    else
        allResult{2}=-9999; %No reflectance grid summary
    end
    allResult{3}=[dorsal_scale,ventral_scale]; %scale length of dorsal and ventral images
    allResult{4}={[both_sides_morph{1}{2},'_',vdlist{both_sides_morph{1}{3}},both_sides_morph{1}{4},'_morph-seg.mat'],...
    [both_sides_morph{2}{2},'_',vdlist{both_sides_morph{2}{3}},both_sides_morph{2}{4},'_morph-seg.mat'],matInNames{1},matInNames{2}}; %original file info
    allResult{5}=tail_morph; %All information related to the tails
    allResult{6}=body_ant_morph; %Copied information related to body size and antennae
    if  procedureFlag<2
        allResult{7}={antennaLRref_rescale, antennaLRrefmm, sppBothSidesAntSummary}; %antennae reflectance
    else
        allResult{7}=-9999; %No antennae reflectance
    end
    %Save all results
    matoutname=fullfile(Result_directory,'dorsal_ventral_map',subFolderList{3},[barcode,'_res-',num2str(2^numberOfIntervalDegree),'x',num2str(2^numberOfIntervalDegree),'_d-v_gridsPars',processingFlagList{procedureFlag+1},'.mat']);
    save(matoutname,'allResult'); %save the specimen matrix
    disp(['################################']);
    disp(['################################']);
    disp(['A matrix of [',barcode,'_res-',num2str(2^numberOfIntervalDegree),'x',num2str(2^numberOfIntervalDegree),'] has been saved' ]);
    disp(['################################']);
    disp(['################################']);
    
    if procedureFlag==0 || procedureFlag==2
        finishedDir='done_d-v';
        if ~exist(fullfile(morph_mat_directory,finishedDir), 'dir')
            mkdir(fullfile(morph_mat_directory,finishedDir));
        end
        movefile(fullfile(morph_mat_directory,[barcodein,'*morph-seg.mat']),fullfile(morph_mat_directory,finishedDir));
    end
catch
    disp('Begin to save all summary matrices WITHOUT reflectance grids or antennae reflectance.');
    allResult=cell(0,4);
    allResult{1}={gridsParameter_dorsal,gridsParameter_ventral}; %All parameters and grids
    allResult{2}=-9999; %No reflectance grid summary
    allResult{3}=[dorsal_scale,ventral_scale]; %scale length of dorsal and ventral images
    allResult{4}={[both_sides_morph{1}{2},'_',vdlist{both_sides_morph{1}{3}},both_sides_morph{1}{4},'_morph-seg.mat'],...
    [both_sides_morph{2}{2},'_',vdlist{both_sides_morph{2}{3}},both_sides_morph{2}{4},'_morph-seg.mat'],'No_reflectance-band_matrix','No_reflectance-band_matrix'}; %original file info
    allResult{5}=tail_morph; %All information related to the tails
    allResult{6}=body_ant_morph; %Copied information related to body size and antennae
    allResult{7}=-9999; %No antennae reflectance
    %Save all results
    matoutname=fullfile(Result_directory,'dorsal_ventral_map',subFolderList{3},[barcode,'_res-',num2str(2^numberOfIntervalDegree),'x',num2str(2^numberOfIntervalDegree),'_d-v_gridsPars-no_reflectance.mat']);
    save(matoutname,'allResult'); %save the specimen matrix
    disp(['################################']);
    disp(['################################']);
    disp(['A matrix of [',barcode,'_res-',num2str(2^numberOfIntervalDegree),'x',num2str(2^numberOfIntervalDegree),'] WITHOUT reflectance data has been saved' ]);
    disp(['################################']);
    disp(['################################']);
   
end

    %Ouput data sturcture:
      %{1}=All parameters and grids
            %{1}{1}=dorsal
            %{1}{2}=ventral
                % {1}{}{1}={seg4PtsLF,wingGridsLF }; %4 key points & grids
                % {1}{}{2}={seg4PtsRF,wingGridsRF}; %4 key points & grids
                % {1}{}{3}={outSeg4PtsLH ,wingGridsLH}; %4 key points & grids
                % {1}{}{4}={outSeg4PtsRH ,wingGridsRH}; %4 key points & grids
                % {1}{}{5}={refineAreaLH,regPtLH,reconstructRegPtLH}; %mask without ornament, 3 key points & key points of the mask without ornament
                % {1}{}{6}={refineAreaRH,regPtRH,reconstructRegPtRH}; %mask without ornament, 3 key points & key points of the mask without ornament
                % {1}{}{7}=original mask of LH
                % {1}{}{8}=original mask of RH
     %{2}=All reflectance grids summary; if it's not applicable, use -9999 as a placeholder
            %{2}{1}=dorsal
            %{2}{2}=ventral
                %{2}{}{1}=740
                %{2}{}{2}=940
                %{2}{}{3}=UV
                %{2}{}{4}=UVF
                %{2}{}{5}=F
                %{2}{}{6}=white
                %{2}{}{7}=whitePo1
                %{2}{}{8}=whitePo2
                %{2}{}{9}=FinRGB
                %{2}{}{10}=PolDiff
                    %{2}{}{}{1}=left forewing
                    %{2}{}{}{2}=right forewing
                    %{2}{}{}{3}=left hindwing
                    %{2}{}{}{4}=right hindwing
                        %{2}{}{}{}{1}=mean reflectance in a grid (gray scale: single layer; RGB: 3 layers)
                        %{2}{}{}{}{2}=standard error of reflectance in a grid (gray scale: single layer; RGB: 3 layers)
      %{3}=scale length of dorsal and ventral images
            %{3}(1)=dorsal
            %{3}(2)=ventral
      %{4}=original file information
            %{4}(1)=dorsal side, morph-seg file name
            %{4}(2)=ventral side, morph-seg file name
            %{4}(3)=dorsal side, AllBandsMask file name
            %{4}(4)=ventral side, AllBandsMask file name
        %{5}=tails information
            %{5}{1}=dorsal_LH tails (-9999 represents no tail)
            %{5}{2}=ventral_LH tails (-9999 represents no tail)
            %{5}{3}=dorsal_RH tails (-9999 represents no tail)
            %{5}{4}=ventral_RH tails (-9999 represents no tail)
                %{5}{}{N}= No. N tail part
                    %{5}{}{}{1}= tail mask
                    %{5}{}{}{2}= raw coordinations of tail base (the part connect to hindwing; points of two ends and the mid-point are provided [mid-point is at second row])
                    %{5}{}{}{3}= summary grid coordinations of tail base (the part connect to hindwing; points of two ends and the mid-point are provided [mid-point is at second row])
                    %{5}{}{}{4}= [tail-base boundary Length, tail length, tail area, tail width, tailCurvature] (in cm; curvature is unit less)
        %{6}= information related to body size and antennae
            %{6}{1}=body length and width (in cm) (dorsal and ventral sides are in two rows)
            %{6}{2}= Antenna length, width, bulb width, degree of curvature (all in mm); LEFT antenna; first row is dorsal side, second is ventral side.
            %{6}{3}= Antenna length, width, bulb width, degree of curvature (all in mm); RIGHT antenna; first row is dorsal side, second is ventral side.
        %{7}= All reflectance on antennae; if it's not applicable, use -9999 as a placeholder
            %{7}{1}= Mean reflectance rescale to 100 grids from the tip of antenna to its base
            %{7}{2}= Mean reflectance every 0.1 mm from the tip of antenna to its base
                %{7}{1|2}{1}= Left antenna
                %{7}{1|2}{2}= Right antenna
                    %{7}{1|2}{}(: , : , 1)= dorsal side
                    %{7}{1|2}{}(: , : , 2)= ventral side
                    %{7}{1|2}{}(: , 1 , :)= distance to the tip of antenna (mm)
                    %{7}{1|2}{}(: , 2:21, :)= spectral reflectance (740, 940, UV, UVF, F, white (R, G, B), whitePo1 (R, G, B), whitePo2 (R, G, B), FinRGB (R, G, B), PolDiff (R, G, B))
            %{7}{3}= Raw reflectance value (raw pixel-based value)
                %{7}{3}{1}= Left antenna
                %{7}{3}{2}= Right antenna
                    %{7}{3}{}{1}= reflectance matrix
                        %{7}{3}{}{1}(:,1)= Percentage from the tip
                        %{7}{3}{}{1}(:,2:21)= spectral reflectance (740, 940, UV, UVF, F, white (R, G, B), whitePo1 (R, G, B), whitePo2 (R, G, B), FinRGB (R, G, B), PolDiff (R, G, B))
                        %{7}{3}{}{1}(:,22)= pixel distance from the tip
                        %{7}{3}{}{1}(:,23)= original coordination X
                        %{7}{3}{}{1}(:,24)= original coordination Y
                    %{7}{3}{}{2}= scale (= 1 cm)
end