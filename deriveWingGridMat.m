function comMatrix=deriveWingGridMat(gapPixel,gridSummaryIn_LF,gridSummaryIn_RF,gridSummaryIn_LH,gridSummaryIn_RH)
%     gridSummaryIn_LF=gridSummaryMeanSELF{1};
%     gridSummaryIn_RF=gridSummaryMeanSERF{1};
%     gridSummaryIn_LH=gridSummaryMeanSELH{1};
%     gridSummaryIn_RH=gridSummaryMeanSERH{1};        

    matrixLF=gridSummaryIn_LF;
    matrixRF=gridSummaryIn_RF;
    matrixLH=gridSummaryIn_LH;
    matrixRH=gridSummaryIn_RH;

%     gapPixel=5;
    comMatrix=vertcat(horzcat(matrixLF, zeros(size(matrixLF,1),gapPixel,size(matrixLF,3)),matrixRF),...
    zeros(gapPixel,2*size(matrixLF,1)+gapPixel,size(matrixLF,3)),...
    horzcat(matrixLH, zeros(size(matrixLH,1),gapPixel,size(matrixLH,3)),matrixRH));
end