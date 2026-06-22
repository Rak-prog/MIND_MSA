function [coordsRAS, labelsOut] = define_schaefer_coord(roi_labels, NPARC, NNET)

name = sprintf("Schaefer2018_%dParcels_%dNetworks_order_FSLMNI152_1mm.Centroid_RAS.csv", NPARC, NNET);
csvPath = fullfile("/home/riccardo/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/MNI/Centroid_coordinates", name);
schaefer_mni = readtable(csvPath);

list1 = string(roi_labels(:));
list1 = strrep(list1, "lh_", "");
list1 = strrep(list1, "rh_", "");

[tf, idx] = ismember(list1, string(schaefer_mni.ROIName));
if any(~tf)
    missing = list1(~tf);
    error("Some labels not found in Schaefer centroid CSV. Example: %s", missing(1));
end

R = schaefer_mni.R(idx);
A = schaefer_mni.A(idx);
S = schaefer_mni.S(idx);

coordsRAS = [R, A, S];
labelsOut = list1; % cleaned ROIName order aligned to roi_labels

end
