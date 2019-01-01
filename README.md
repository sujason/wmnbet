# wmnbet: White-matter-nulled (MPRAGE) brain extraction tool
Segmentation of the brain parenchyma and intra-cranial cavity (ICC) using the white-matter nulled contrast and label fusion.  Note this requires label_fusion_data/, provided elsewhere.

This tool was created based on data acquired on ET FUS subjects.  The template is formed from 3T WMnMPRAGE images.  Prior brain parenchyma segmentations were created by applying BET to CSFnMPRAGE images co-registered to WMnMPRAGE.  Prior ICC segmentations were created by applying BET to T2CUBE images co-registered to WMnMPRAGE.  These sequences were acquired in the same scan session and co-registered using linear registration.

## Requirements
- [ANTs](https://github.com/stnava/ANTs.git)
	- NOTE: before April 2016, there was a bug in ANTs/Examples/antsJointFusion.cxx when using the "-x" switch, please rebuild ANTs from the latest source to fix this.
- [FSL](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation)
- [convert3d](http://www.itksnap.org/pmwiki/pmwiki.php?n=Downloads.C3D)

## Usage
- python label_fusion.py check brain.json
- python label_fusion.py do --tempdir temp --labelfusion average -p 2 brain.json image.nii.gz output/ brain icc
