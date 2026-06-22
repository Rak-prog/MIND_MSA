function define_DK_coord_and_save(roi_labels, output_path)

    dk_mni = readtable('/mnt/NAS_Progetti/PSP-MSA/daPellecchia/Bids_Salerno/derivatives/structural_similarity/network_statistics_NBS_old/ROI_coord_DK.txt');
    dk_mni.labelname = strcat(dk_mni.Hemisphere, '_', dk_mni.ParcelName);
    dk_mni.labelname = strrep(dk_mni.labelname , 'left', 'lh');
    dk_mni.labelname = strrep(dk_mni.labelname , 'right', 'rh');
    
    dk_mni.ParcelName = [];
    dk_mni.Hemisphere = [];
    
    
    % fuzzy matching string with coordinate and the one caluclated with MIND to
    % get the MNI roi needed for Network Based Statistics 
    % devo ordinare le coordinate MNI a roi_labels che sono le matrici di
    % adiacenza, IMPORTANTE :) 
    list1 = roi_labels';  % this are the coordinates used in the the MIND code, so I align the coordinates to this one
    [~, idx] = ismember(list1, dk_mni.labelname);
    % Handle missing labels (if any)
    if any(idx == 0)
        list1(idx==0)
        error('Some labels in list1 were not found in dk_mni.labelname');
    end
    % Reorder dk_mni according to list1
    dk_mni_reordered = dk_mni(idx, :);
    check_table = table(list1(:), dk_mni_reordered.labelname(:), ...
        'VariableNames', {'DesiredOrder', 'MatchedLabel'});
    % Save results
    writematrix(string(dk_mni_reordered.labelname), fullfile(output_path, 'ROI_labels_final.txt'));
    writematrix([dk_mni_reordered.MNI_X, dk_mni_reordered.MNI_Y, dk_mni_reordered.MNI_Z], ...
                 fullfile(output_path, 'ROI_coordMNI_final.txt'));

end