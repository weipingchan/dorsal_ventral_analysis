function saveTemporaryUniqueBarcodeListOfMorphSeg(morph_mat_directory,Code_directory)
    addpath(genpath(Code_directory)) %Add the library to the path
%     all_morph_data = dir(fullfile(morph_mat_directory,['**',filesep,'*morph-seg.mat']));
    all_morph_data = dir(fullfile(morph_mat_directory,'*morph-seg.mat'));
    %generate a unique barcode list
    barcodelist=cell(0,1);
    rec=1;
    for matinID=1:length(all_morph_data)
        matinname=all_morph_data(matinID).name;
        [barcode, side, flag]=file_name_decoder(matinname);

        %         if  isempty(find(contains(barcodelist,barcode)))
        if isempty(find(strcmp(barcodelist,barcode)))
            barcodelist{rec}=barcode;
            rec=rec+1;
        end
    end

    fid = fopen(fullfile(morph_mat_directory,'tmp_unique_barcode_list.txt'), 'w');
    for row = 1:length(barcodelist)
        fprintf(fid, '%s\n', barcodelist{row});
    end
    fclose(fid);
end