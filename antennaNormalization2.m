function [antennaLRref_rescale, antennaLRrefmm]=antennaNormalization2(sppBothSidesAntSummary,outPutLevel)
% outPutLevel=100;
    antennaLRrefmm=cell(0,2);
    antennaLRref_rescale=cell(0,2);
    if ~isempty(sppBothSidesAntSummary)
        for LRid=1: 2
            if ~isempty(sppBothSidesAntSummary{LRid})
                antennae_d=sppBothSidesAntSummary{LRid}{1};
                antennae_v=sppBothSidesAntSummary{LRid}{2};

                [antMatmm_d, emptyFlag2_d] = normalize_antMat(antennae_d);
                [antMatmm_v, emptyFlag2_v] = normalize_antMat(antennae_v);
                sizeCand=[size(antMatmm_d,1),size(antMatmm_v,1)];
                sizeCand=sizeCand(sizeCand>20); %Antennae should be at least 2 mm long
                keepLen=min(sizeCand); %Keep data from the tips of antennae in both dorsal and ventral sides
                
                if emptyFlag2_d==1
                    antMatmm_d=zeros(keepLen,21);
                    antMatmm_d(antMatmm_d==0)=NaN;
                end
                if emptyFlag2_v==1
                    antMatmm_v=zeros(keepLen,21);
                    antMatmm_v(antMatmm_v==0)=NaN;
                end
                
                if emptyFlag2_d~=1
                    antMatmm_d2=antMatmm_d(1:keepLen,:);
                    rescale_antmm_d2=imresize(antMatmm_d2,[outPutLevel,size(antMatmm_d2,2)]);                
                else
                    rescale_antmm_d2=zeros(outPutLevel,size(antMatmm_d,2));
                    rescale_antmm_d2(rescale_antmm_d2==0)=NaN;
                end
                if emptyFlag2_v~=1
                    antMatmm_v2=antMatmm_v(1:keepLen,:);
                    rescale_antmm_v2=imresize(antMatmm_v2,[outPutLevel,size(antMatmm_v2,2)]);
                else
                    rescale_antmm_v2=zeros(outPutLevel,size(antMatmm_v,2));
                    rescale_antmm_v2(rescale_antmm_v2==0)=NaN;
                end
                antMatmm_dv=cat(3, antMatmm_d(1:keepLen,:),antMatmm_v(1:keepLen,:));
                antMatmm_dv_rescale=cat(3, rescale_antmm_d2,rescale_antmm_v2);
            else
                antMatmm_dv=[];
                antMatmm_dv_rescale=[];
            end
            antennaLRrefmm{LRid}=antMatmm_dv;
            antennaLRref_rescale{LRid}=antMatmm_dv_rescale;
        end
    end

end