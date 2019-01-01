# wmnbet: White-matter-nulled (MPRAGE) brain extraction tool
Segmentation of the brain parenchyma and intra-cranial cavity (ICC) using the white-matter nulled contrast and PICSL's joint label fusion.  Note this requires label_fusion_data/, provided elsewhere.

This tool was created based on data acquired on ET FUS subjects.  The template is formed from 3T WMnMPRAGE images.  Prior brain parenchyma segmentations were created by applying BET to CSFnMPRAGE images co-registered to WMnMPRAGE.  Prior ICC segmentations were created by applying BET to T2CUBE images co-registered to WMnMPRAGE.  These sequences were acquired in the same scan session and co-registered using linear registration.
