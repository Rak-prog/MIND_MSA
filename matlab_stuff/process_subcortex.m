function subcort_final = process_subcortex(subcort)
    subcort_final = table();
    subcort_final.Putamen = subcort.leftPutamen + subcort.rightPutamen;
    subcort_final.CerebWM = subcort.leftCerebellumWhiteMatter + subcort.rightCerebralWhiteMatter;
    subcort_final.CerebGM = subcort.leftCerebellumCortex + subcort.rightCerebellumCortex;
    subcort_final.Pons = subcort.Pons;
    subcort_final.Midbrain = subcort.Midbrain;
    subcort_final.Pallidum = subcort.leftPallidum + subcort.rightPallidum;
    subcort_final.Caudate = subcort.leftCaudate + subcort.rightCaudate;
    subcort_final.Thalamus = subcort.leftThalamus + subcort.rightThalamus;
    subcort_final.Whole_brainstem = subcort.Whole_brainstem;
    subcort_final.ICV = subcort.totalIntracranial;
    subcort_final.Subj_ID =  replace(subcort.Subj_ID, "_ses-T0", "");
end