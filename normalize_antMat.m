function [antMatmm, emptyFlag2] = normalize_antMat(antennae)
    emptyFlag2=0;
    if ~isempty(antennae)
        if ~isempty(antennae{1}) && size(antennae{1},2)>1
            scale=antennae{2}{1};
            antMat=sort(antennae{1});
             [uniqueA i j] = unique(antMat(:,22),'first');
             indexToDupes = find(not(ismember(1:numel(antMat(:,22)),i))); %Find duplicate value
             antMat(indexToDupes,22)=antMat(indexToDupes,22)+0.0001; %Make it slightly different
            antL=max(antMat(:,22))/scale*10; %mm based on the new measurement
            intpLoc = interp1(antMat(:,22)/scale*10 , [1:1:length(antMat(:,22))] , [0 : 0.1 : antL]);
            roundintpLoc = round(intpLoc);
            refintmm=zeros(size(intpLoc,2)-1,20);
            for intID=1:size(intpLoc,2)-1
                refInt=mean(antMat(roundintpLoc(intID):roundintpLoc(intID+1),2:21),1);
                refintmm(intID,:)=refInt;
            end
            antMatmm=[[0.1 : 0.1 : antL]', refintmm];
            if size(antMatmm,1)<=20  %Antennae should be at least 2 mm long
                emptyFlag2=1;
            end
        else
            antMatmm=[];
            emptyFlag2=1;
        end
    else    
        antMatmm=[];
        emptyFlag2=1;
    end
end