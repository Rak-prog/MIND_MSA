function masks = build_class_masks(conn2, demo2)

grp = upper(strtrim(string(conn2.group)));

masks.HC  = grp=="HC";
masks.PD  = grp=="PD";
masks.MSA = grp=="MSA";
masks.PSP = grp=="PSP";
masks.HCyoung = grp=="HCYOUNG";
% --- MSA subtype (edit demo column name) ---
msa_sub = upper(strtrim(string(demo2.Phenotype)));
masks.MSAP = masks.MSA & contains(msa_sub,"P");
masks.MSAC = masks.MSA & contains(msa_sub,"C");
if all(conn2.group == demo2.GROUP) == 1
    masks.MSA_P1 = grp=="MSA" & demo2.PROTOCOL == 1;
    masks.MSA_P2 = grp=="MSA" & demo2.PROTOCOL == 2;
    masks.MSA_P3 = grp=="MSA" & demo2.PROTOCOL == 3;
    masks.MSAC_P1 = masks.MSA & contains(msa_sub,"C") & demo2.PROTOCOL == 1;
    masks.MSAC_P2 = masks.MSA & contains(msa_sub,"C") & demo2.PROTOCOL == 2;
    masks.MSAC_P3 = masks.MSA & contains(msa_sub,"C") & demo2.PROTOCOL == 3;
    masks.MSAP_P1 = masks.MSA & contains(msa_sub,"P") & demo2.PROTOCOL == 1;
    masks.MSAP_P2 = masks.MSA & contains(msa_sub,"P") & demo2.PROTOCOL == 2;
    masks.MSAP_P3 = masks.MSA & contains(msa_sub,"P") & demo2.PROTOCOL == 3;
else
    error("demo2 and conn matrices are not matching")
end


% --- MCI within MSA (edit demo column name) ---
mci = demo2.MCI_T0;
if ~isnumeric(mci)
    mci = upper(strtrim(string(mci)));
    mci = double(ismember(mci, ["YES","Y","1","TRUE"]));
end
masks.MCIyes = masks.MSA & (mci==1);
masks.MCIno  = masks.MSA & (mci==0);

% --- PSP subtypes (edit demo column name) ---
psp_sub = upper(strtrim(string(demo2.Phenotype)));
masks.PSP_cor = masks.PSP & (contains(psp_sub, "PSP-F") | contains(psp_sub, "PSP-CBS") | contains(psp_sub, "PSP-CBD"));
masks.PSP_sub = masks.PSP & (contains(psp_sub, "PSP-P") | contains(psp_sub, "PSP-PGF") | contains(psp_sub, "PSP-OM"));
masks.PSP_R   = masks.PSP & contains(psp_sub, "PSP-RS");

end
