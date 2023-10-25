# BIDS Aging Data
**GOAL:** 

BIDSification of aging paper's processed data ([Callaghan et al., 2014](https://doi.org/10.1016/j.neurobiolaging.2014.02.008)) for public release.

Some aspects of the data set description are discussed in this [HackMd note](https://hackmd.io/@cphillips/B1jhtqVCn) and here we only present

- some resources and arguments
- the resulting BIDS formatting proposition

---

## Data & information available

### qMRI data

Data used for the VBQ analysis
- 138 subjects, age range [min max], 49 M and 89 F
- for each subject, 8 qMR images in normalized space
    - GM-smoothed A (like PD), MT, R1 and R2* maps
    - WM-smoothed A (like PD), MT, R1 and R2* maps
- mutually exclusive binary masks for GM and WM tissues
- 4 regressors with values for subject's age, gender, TIV & scanner 
    - age is expressed in years;
    - gender is 1 for 'male', 0 for 'female';
    - TIV (total intracranial volume) is expressed in liters;
    - scanner is 0 for 'quatro' and 1 for 'trio' but also encoded in the filename with `MQ` (for 'quatro') or `MT` (for 'trio'), with 3 subjects using lower case labels `mq`/`mt`.

Data used for the VBM analysis
- smoothed modulated normalized GM maps from the 138 subjects
- same 4 regressors as for VBQ
- maybe a brain mask?

Then for display purpose a mean MT image but other mean images (R1) are possible.

***Questions***
- What about the smooth modulated warp GM density maps, used for the VBM analysis? :arrow_forward: it would be nice for the sake of completeness.
- Should be included as they are part of the original paper? :arrow_forward: CP: we do not have these at the CRC.
- Should we randomize the subjects order and ensuing BIDS-labels? :arrow_forward: yes, it does not cost much in doing so and improves anonymization (though not perfect...)

### Previous preprocessing steps

Main processing steps that took place to produce the shared preprocessed data
1. creation of the quantitative maps with the "VBQ toolbox", now "hMRI toolbox" :arrow_forward: A, MT, R1 and R2* maps
2. segmentation of MT maps :arrow_forward: GM & WM density maps, in subject space + "DARTEL imported"
3. estimation, based on individual GM and WM maps, of the warping into a MNI-like population specific template space
4. warping of the quantitative (A, MT, R1 and R2*) maps, no modulation
5. warping of the GM & WM density maps, with modulation
6. Gaussian smoothing of the warped modulated GM density maps :arrow_forward: data for VBM analysis
7. GM and WM tissue-weighted smoothing of the warped quantitative maps :arrow_forward: data for VBQ analysis
8. creation of tissue specific GM & WM mask using a "winner-takes-it-all" approach from the mean of all smooth warped modulated GM & WM (+CSF?) density maps :arrow_forward: masks for VBQ analysis
9. creation of mean warped MT maps :arrow_forward: background for VBQ results display

***Questions***
- These should be at least briefly described in the meta-data, as provenance stuff? 
- was CSF used in the GM/WM masks creation? it should have.
- Just assuming data were processed the "standard way" :arrow_forward: to be checked

### References & resources

Here we'll need to improvise keeping BIDS principle in mind... Here are some resources:
- [BIDS starter kit](https://bids-standard.github.io/bids-starter-kit/folders_and_files/derivatives.html)
- [BIDS manual for derivatives](https://bids-specification.readthedocs.io/en/stable/05-derivatives/01-introduction.html)
- [BEP11 for structural preprocessing derivatives](https://docs.google.com/document/d/1YG2g4UkEio4t_STIBOqYOwneLEs1emHIXbGKynx7V0Y/edit?usp=sharing) :arrow_forward: useful for the GM & WM maps ?
- [BIDS for quantitative MRI](https://bids-specification.readthedocs.io/en/stable/appendices/qmri.html) :arrow_forward: useful for the quantitative maps
- [BIDS Extension Proposal 38 (BEP038): Atlas Specification](https://docs.google.com/document/d/1RxW4cARr3-EiBEcXjLpSIVidvnUSHE7yJCUY91i5TfM/edit#heading=h.4k1noo90gelw) :arrow_forward: useful for the tissue masks?

---

## Proposed naming scheme for aging dataset

Overall there are (at least) 3 types of images that need to be described: 
- smoothed warped tissue density maps
- tissue masks & average maps
- TW-smoothed warped quantitative maps

Given their different nature, the naming suffixes will be quite different.

Bringing everything together, we could end up with the following naming scheme:
- at the "study" level, i.e. whole derivatives data, some "atlases" for the masks and mean images:
````
AgingData/derivatives/SPM12dartel/
   atlas-GM_space-MNI_mask.nii
   atlas-GM_space-MNI_mask.json
   atlas-WM_space-MNI_mask.nii
   atlas-WM_space-MNI_mask.json
   atlas-MTsat_space-MNI_desc-mean.nii
   atlas-MTsat_space-MNI_desc-mean.json
   atlas-PDmap_space-MNI_desc-mean.nii
   atlas-PDmap_space-MNI_desc-mean.json
   atlas-R1map_space-MNI_desc-mean.nii
   atlas-R1map_space-MNI_desc-mean.json
   atlas-R2starmap_space-MNI_desc-mean.nii
   atlas-R2starmap_space-MNI_desc-mean.json

````
- at the "study" level, i.e. whole derivatives data, some extra`.json` files about describing the tissue-weighted smoothed warped images from each subjects:
````
AgingData/derivatives/SPM12dartel/
   MTsat_space-MNI_desc-TWsmo.json
   PDmap_space-MNI_desc-TWsmo.json
   R1map_space-MNI_desc-TWsmo.json
   R2starmap_space-MNI_desc-TWsmo.json
````
- at the "subject" level, i.e. in each subject's `anat` subfolder, all the smooth  warped individual images:
````
AgingData/derivatives/SPM12dartel/
   sub-S123/anat/
      sub-S123_MTsat_space-MNI_desc-smomod_label-GM_probseg.nii
      sub-S123_MTsat_space-MNI_desc-smomod_label-WM_probseg.nii
      sub-S123_MTsat_space-MNI_desc-smomod_probseg.json
      sub-S123_MTsat_space-MNI_desc-GMsmo.nii
      sub-S123_MTsat_space-MNI_desc-WMsmo.nii
      sub-S123_PDmap_space-MNI_desc-GMsmo.nii
      sub-S123_PDmap_space-MNI_desc-WMsmo.nii
      sub-S123_R1map_space-MNI_desc-GMsmo.nii
      sub-S123_R1map_space-MNI_desc-WMsmo.nii
      sub-S123_R2starmap_space-MNI_desc-GMsmo.nii
      sub-S123_R2starmap_space-MNI_desc-WMsmo.nii
````

---

## Original processed data organization 

The original processed data  are organized like this, on my HD from the top folder with path `C:\Dox\2_Data\qMRI_MPM\Data4ChrisPhilips`:

- a `.mat` file (`Subjects4Chris.mat`) with ID & regressors (age, TIV, scanner, gender);
- two mask files, `WinnerTakesAllMask_GM.nii` and `WinnerTakesAllMask_WM.nii`, for the GM and WM, resolution `[1 1 1]`mm^3;
- the `Average_MTsat_Maps` folder with 3 average MT images
  - `Averaged_MT_HeadMask.nii`, average across all subjects of whole warped MTmap images, resolution `[1 1 1]`mm^3,  masked including some non-brain tissues;
  - `AveragedMT_MNI_BrainMask.nii`, average across all subjects of whole warped MTmap images, resolution `[1 1 1]`mm^3,  masked tight around the brain (but not matching the GM/WM masks!);
  - `AveragedMTMap_DARTEL.nii`, average across all subjects of whole warped MTmap images, resolution `[1.5 1.5 1.5]`mm^3, no masking;
- the `Fin_dart_p1` folder, with GM-weighted warped maps used for the statistical analysis
  - 3 files `fin_dart_p1MQ0190_R1.nii`, `fin_dart_p1MQ0190_R2s.nii` &  `fin_dart_p1MQ0190_A.nii` with the averaged across all subjects of the `R1`, `R2s` & `A` GM-weighted warped maps;
  - 4 folders `Imgs_A`, `Imgs_MT`, `Imgs_R1 ` & `Imgs_R2s` with all the subjects `A`, `MT` `R1` & `R2s` GM-weighted warped maps
  - each of the 4 `A`, `MT` `R1` & `R2s` folders contains files of the form `fin_dart_p1MT02253_MT.nii` where
    - `MT02253` encodes for the subjects index, `02253`, and scanner label, `MT` (other scanner would be `MQ`);
    - `_MT` for MTsat images, there are thus corresponding images with the `A`, `R1` & `R2s` suffixes.
- the same for the `Fin_dart_p2` folder, with WM-weighted warped maps used for the statistical analysis










````
dsfg
````

