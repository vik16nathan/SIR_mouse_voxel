#!/bin/bash
# Required: minc-toolkit, ANTS, pyminc
# Usage:
# align_to_Allen.sh [Your MR template to align to CCFv3] [Your MR template mask] [Output base name (file extensions will be automatically added)] [Output directory] [Resolution of Allen template, pick one of: 10, 25, 50, 100]

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREAD=8


INFILE=$1 #Image to align to template
INMASK=$2 #Input image brain mask
OUTBASE=$3 # Output basename
OUTDIR=$4 #Directory to which to output the data
RESOLUTION=$5 # Resolution of the Allen template to use

# File and directory locations
# Can change
TFM_MICE=/data/chamal/projects/yohan/common/allenbrain/code/scripts/transform_space.py

AFFINEDIR=${OUTDIR}/affine
NLINDIR=${OUTDIR}/nlin
TEMPLATE=${OUTDIR}/average_template.mnc
TEMPLATE_EXPANDED=${OUTDIR}/average_template_expanded.mnc
TEMPLATE_MASK=${OUTDIR}/average_template_mask.mnc
TEMPLATE_MASK_EXPANDED=${OUTDIR}/average_template_mask_expanded.mnc
INFMASKED=${OUTDIR}/input_masked.mnc

AFFINE=${AFFINEDIR}/${OUTBASE}_affine.mnc
AFFINE_MASK=${AFFINEDIR}/${OUTBASE}_affine_mask.mnc
AFFINE_MASKED=${AFFINEDIR}/${OUTBASE}_affine_masked.mnc
NLIN=$NLINDIR/${OUTBASE}_nlin.mnc
NLIN_MASK=${NLINDIR}/${OUTBASE}_nlin_mask.mnc
NLIN_MASKED=${NLINDIR}/${OUTBASE}_nlin_masked.mnc


#URL and file definitions
# Don't change
DOWNLOAD_TEMPLATE="average_template_${RESOLUTION}.nrrd"
DOWNLOAD_MASK="structure_997.nrrd"
TEMPLATE_URL="https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/average_template/${DOWNLOAD_TEMPLATE}"
TEMPLATE_MASK_URL="https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/structure_masks/structure_masks_${RESOLUTION}/${DOWNLOAD_MASK}"

echo "Getting input files..."
mkdir ${OUTDIR}
cp ${INFILE} ${OUTDIR}/input_image.mnc
cp ${INMASK} ${OUTDIR}/input_mask.mnc
wget -P ${OUTDIR} ${TEMPLATE_URL} --no-check-certificate
wget -P ${OUTDIR} ${TEMPLATE_MASK_URL} --no-check-certificate

# Convert to MINC
echo "Converting files..."
eval "${TFM_MICE} ${OUTDIR}/${DOWNLOAD_TEMPLATE} ${TEMPLATE} -v RAS -w MICe -x 1"
eval "${TFM_MICE} ${OUTDIR}/${DOWNLOAD_MASK} ${TEMPLATE_MASK} -v RAS -w MICe -x 1"
if [ "${RESOLUTION}" -eq "100" ]; then
    mincmorph -successive DDDDDDEEEEEEDD ${TEMPLATE_MASK} ${TEMPLATE_MASK_EXPANDED}
fi
if [ "${RESOLUTION}" -eq "50" ]; then
    mincmorph -successive DDDDDDDDDDDDEEEEEEEEEEEEDDDD ${TEMPLATE_MASK} ${TEMPLATE_MASK_EXPANDED}
fi
if [ "${RESOLUTION}" -eq "25" ]; then
    mincmorph -successive DDDDDDDDDDDDDDDDDDDDDDDDEEEEEEEEEEEEEEEEEEEEEEEEDDDDDDDD ${TEMPLATE_MASK} ${TEMPLATE_MASK_EXPANDED}
fi

mincmask ${INFILE} ${INMASK} ${INFMASKED}
rm ${OUTDIR}/${DOWNLOAD_TEMPLATE} ${OUTDIR}/${DOWNLOAD_MASK}

echo "Expanding template..."
autocrop -isoexpand 20% ${TEMPLATE} ${TEMPLATE_EXPANDED}

# Affine registration
echo "Affine registration..."
mkdir ${AFFINEDIR}

CMD="/usr/bin/time -v antsRegistration \
-d 3 \
-o ${AFFINEDIR}/${OUTBASE}_affine \
-a 1 \
-z 1 \
--initial-moving-transform [${INFILE},${TEMPLATE},1] \
-t Affine[0.25] \
-m MI[${INFILE},${TEMPLATE},1,32,Regular,0.25] \
--convergence [1000x500x250x100,1e-6,10] \
--shrink-factors 8x4x2x1 \
--smoothing-sigmas 0.2x0.1x0.05x0mm \
--minc \
-v 1 \
--use-histogram-matching 1"

echo "$CMD" > ${OUTDIR}/ANTS_affine.cmd
eval ${CMD}

mincresample -2 -like ${TEMPLATE_EXPANDED} -transform ${AFFINEDIR}/${OUTBASE}_affine.xfm ${INFILE} ${AFFINE}
mincresample -2 -like ${TEMPLATE_EXPANDED} -transform ${AFFINEDIR}/${OUTBASE}_affine.xfm ${INFMASKED} ${AFFINE_MASKED}
mincresample -2 -like ${TEMPLATE_EXPANDED} -transform ${AFFINEDIR}/${OUTBASE}_affine.xfm -nearest_neighbour -keep_real_range ${INMASK} ${AFFINE_MASK}

# Nonlinear registration
echo "Nonlinear registration"
mkdir ${NLINDIR}

CMD="/usr/bin/time -v antsRegistration \
-d 3 \
-o ${NLINDIR}/${OUTBASE}_nlin \
-a 1 \
-z 1 \
-t SyN[0.1,2,0] \
-m MI[${AFFINE},${TEMPLATE},1,32,Regular,0.25] \
--convergence [1000x800x500x200,1e-6,10] \
--shrink-factors 8x4x2x1 \
--smoothing-sigmas 0.2x0.1x0.05x0mm \
--minc \
-v 1 \
-x ${TEMPLATE_MASK_EXPANDED} \
--use-histogram-matching 1"

echo "$CMD" > ${OUTDIR}/ANTS_nlin.cmd
eval $CMD

mincresample -2 -like ${TEMPLATE} -transform ${NLINDIR}/${OUTBASE}_nlin.xfm ${AFFINE} ${NLIN}
mincresample -2 -like ${TEMPLATE} -transform ${NLINDIR}/${OUTBASE}_nlin.xfm ${AFFINE_MASKED} ${NLIN_MASKED}
mincresample -2 -like ${TEMPLATE} -transform ${NLINDIR}/${OUTBASE}_nlin.xfm -nearest_neighbour -keep_real_range ${AFFINE_MASK} ${NLIN_MASK}

ln -s ${NLIN} ${OUTDIR}/input_image_resampled.mnc
ln -s ${NLIN_MASK} ${OUTDIR}/input_mask_resampled.mnc
xfmconcat ${AFFINEDIR}/${OUTBASE}_affine.xfm ${NLINDIR}/${OUTBASE}_nlin.xfm ${OUTDIR}/${OUTBASE}.xfm
xfminvert ${OUTDIR}/${OUTBASE}.xfm ${OUTDIR}/${OUTBASE}_inverted.xfm
