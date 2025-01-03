#!/bin/bash
#
#      0        1        2   
#  vizz.sh  <subj_id>  <hemi>
#
#
#

subj_id=$1
hemi=$2

Doreti_dir=/mnt/neurocube/local/serenceslab/Doreti
subj_dir=$Doreti_dir/FUNC/$subj_id/PreProc/AFNI
SUMA_dir=$Doreti_dir/ANAT/$subj_id/SUMA
pal_dir=$Doreti_dir/SCRIPTS_DOTS

portnumber=$(afni -available_npb_quiet) # prints & saves first available port number to a variable
#AFNI_COLORSCALE_DEFAULT=HSV_${hemi}

#overlayfile=params_coords_sm5.nii.gz
#overlayfile=sfaout_polar+orig.BRIK
# overlayfile=tstat_bowtie_H+orig



# suma -niml -npb $portnumber -spec ${SUMA_dir}/${subj_id}_${hemi}.spec -sv ${SUMA_dir}/${subj_id}_SurfVol.nii &
# afni -niml -npb $portnumber -yesplugouts -posfunc -skip_afnirc \
# 	$subj_dir/$overlayfile &
# plugout_drive -npb $portnumber  \
# 	-com 'SET_THRESHOLD A.000'                               \
# 	-com 'SET_FUNC_RANGE A.8'                               \
# 	-com 'SEE_OVERLAY +'       


# overlayfile=sfaout_polar+orig
overlayfile=polar_phz_adj_MMH+orig


suma -niml -npb $portnumber -spec ${SUMA_dir}/${subj_id}_${hemi}.spec -sv ${SUMA_dir}/${subj_id}_SurfVol.nii &
afni -niml -npb $portnumber -yesplugouts -posfunc -skip_afnirc \
	-DAFNI_COLORSCALE_38=$pal_dir/HSV_lh.pal \
	-DAFNI_COLORSCALE_39=$pal_dir/HSV_rh.pal \
	-DAFNI_COLORSCALE_40=$pal_dir/Jet256.pal \
	-DAFNI_COLORSCALE_DEFAULT=HSV_${hemi} \
	$subj_dir/$overlayfile &
plugout_drive -npb $portnumber  \
	-com 'SET_THRESHOLD A.0'                               \
	-com 'SET_FUNC_RANGE A.360'                               \
	-com 'SEE_OVERLAY +'                           

# suma -niml -npb $portnumber -spec ${SUMA_dir}/${subj_id}_${hemi}.spec -sv ${SUMA_dir}/${subj_id}_SurfVol.nii &
# afni -niml -npb $portnumber -yesplugouts -posfunc -skip_afnirc \
# 	-DAFNI_COLORSCALE_38=$pal_dir/HSV_lh.pal \
# 	-DAFNI_COLORSCALE_39=$pal_dir/HSV_rh.pal \
# 	-DAFNI_COLORSCALE_40=$pal_dir/Jet256.pal \
# 	-DAFNI_COLORSCALE_DEFAULT=HSV_${hemi} \
# 	$subj_dir/$overlayfile &
# plugout_drive -npb $portnumber  \
# 	-com 'SET_SUBBRICKS A 0 0 2'                              \
# 	-com 'SET_THRESHOLD A.0100'                               \
# 	-com 'SET_FUNC_RANGE A.360'                               \
# 	-com 'SEE_OVERLAY +'                               
