function sppBothSidesAntSummary=extractAntennaReflectance(spp_mat_directory,both_sides_morph, ants_d, ants_v, dorsal_key, ventral_key_flip, dorsal_scale, ventral_scale)
 sppBothSidesAntSummary=cell(0,2);% Antennae
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
            disp([matinname,' has been read into system and ready for process.']);

             %%
            %Antennae extraction
            disp('Begin to summarize multi-spectral reflectance on antenna');
            mask0=both_sides_morph{side}{1}{1};
            wingParts0=both_sides_morph{side}{1}{7};
            wingParts=wingParts0{1}+wingParts0{2}+wingParts0{3}+wingParts0{4};
            if side==1
                mask=mask0;
                wingMask=wingParts;
                ants_mask=ants_d;
                realCen=dorsal_key(8,:);
                spScale=dorsal_scale;
            elseif side==2
                mask=flip(mask0,2);
                wingMask=flip(wingParts,2);
                ants_mask=ants_v;
                realCen=ventral_key_flip(8,:);
                spScale=ventral_scale;
            end
            

            antenna0=bwareafilt((mask-wingMask)>0,1);
            BWant0 = bwskel(logical(antenna0),'MinBranchLength',30);

            [labeledAntenna, numberOfBlobs] = bwlabel(ants_mask);
            sppBothSidesAntSummary{1}{sideID}={-9999,{spScale}};
            sppBothSidesAntSummary{2}{sideID}={-9999,{spScale}};
            if numberOfBlobs>0
                try
                    if numberOfBlobs>1  %Get antenae Bases
                        antennaeBaseList=[[NaN NaN];[NaN NaN]];
                        antennaeBaseList=[];
                        for k=1:numberOfBlobs
                            try
                            oneAntenae=ismember(labeledAntenna, k);
                            antenaeEdgeMask=imdilate(oneAntenae,strel('disk',1))-imerode(oneAntenae,strel('disk',1));
                            [ei,ej] = find(immultiply(BWant0,antenaeEdgeMask));
                            antennaeBase0=sortrows([ei,ej]); %sort
                            antennaeBase=antennaeBase0(1,:);
                            antennaeBaseList(k,:)=antennaeBase;
                            catch
                            end
                        end
                        if ~isempty(antennaeBaseList)
                            antennaeBaseList=sortrows(antennaeBaseList,2);
                        end
                    end

                    for k=1:numberOfBlobs
                        try
                            oneAntenae=ismember(labeledAntenna, k);
                            antenaeEdgeMask=imdilate(oneAntenae,strel('disk',1))-imerode(oneAntenae,strel('disk',1));
                            [ei,ej] = find(immultiply(BWant0,antenaeEdgeMask));
                            antennaeBase0=sortrows([ei,ej]); %sort
                            antennaeBase1=antennaeBase0(1,:);

                %             figure,imshowpair(BWant0,antenaeEdgeMask);

                            antPath=immultiply(BWant0,oneAntenae);

                            [ti,tj] = find(antPath);

                            antennaeBase=findCloestPt([ti,tj],antennaeBase1);

                            Dd = bwdistgeodesic(antPath,antennaeBase(2),antennaeBase(1),'quasi');
                            coAnt=[];
                            for n = 1:numel(ti)
                                coAnt=[coAnt; [ti(n), tj(n), Dd(ti(n),tj(n))]];
                            end 

                            antAxis=(max(coAnt(:,3))-coAnt(:,3))/max(coAnt(:,3));

                            sppAntRefBands0=[];
                            for bandID=1:10
                                if side==1
                                    baseImg=sppmat{bandID};
                                elseif side==2
                                    baseImg=flip(sppmat{bandID},2);
                                end
                                refAnt=[];
                                for n = 1:numel(ti)
                                    refAnt=[refAnt; [Dd(ti(n),tj(n)), reshape(baseImg(ti(n),tj(n),:),1,[])]];
                                end
                                sppAntRefBands0=[sppAntRefBands0, refAnt(:,2:end)];
                            end

                            sppAntRefBands=[antAxis, sppAntRefBands0, round((max(coAnt(:,3))-coAnt(:,3)), 2), coAnt(:,1:2)];
                            %Percentage from the tip, spectral reflectance, pixel distance from the tip, original coordination (x, y)
                                %spectral reflectance: 740, 940, UV, UVF, F, white (R, G, B), whitePo1 (R, G, B), whitePo2 (R, G, B), FinRGB (R, G, B), PolDiff (R, G, B)

                            %Determine if it's left or right antennae
                            ss01 = regionprops(oneAntenae,'centroid');
                            antCens=ss01.Centroid;
                            if numberOfBlobs<2
                                if antCens(1)<realCen(1)
                                    rfAnt=1;
                                else
                                    rfAnt=2;
                                end
                            else
                                [iloc, ~] = ismember(antennaeBaseList,antennaeBase,'rows'); %If the base of the antennae is at left than it's the left one
                                rfAnt=find(iloc);
                            end

                            if side*rfAnt~=2
                                sppBothSidesAntSummary{1}{sideID}={sppAntRefBands,{spScale}};
                            else
                                sppBothSidesAntSummary{2}{sideID}={sppAntRefBands,{spScale}};
                            end
                %             %Read corresponding parameters of antennae
                %             % L ant            side==1, rfAnt==1
                %             % R ant            side==1, rfAnt==2
                %             % R ant            side==2, rfAnt==1
                %             % L ant            side==2, rfAnt==2
                %             if side*rfAnt~=2
                %                 antChr0=body_ant_morph{2};
                %             else
                %                 antChr0=body_ant_morph{3};
                %             end
                %             antChr=antChr0(side,:); %Antennae length, width, bolb width, degree of curved (all in mm)
                        catch
                            disp('wrong antenna processing');
                        end
                    end
                catch
                    disp('wrong antenna processing');
                end
            else
                    disp('No Antennae');
            %             sppBothSidesAntSummary{1}{sideID}={-9999,{spScale}};
            %             sppBothSidesAntSummary{2}{sideID}={-9999,{spScale}};
            end
            disp([vdlist{sideID},' side of Antennae in [',barcode,'] has been extracted' ]);
        end
end