design=322_fsl_
dir=/vols/Scratch/mgarvert/ManyMaps/imagingData/2ndLevel/design_${design} 
maskdir=/vols/Scratch/mgarvert/ManyMaps/imagingData/masks 
for mask in 322_1_07_left_parahippoc_2p5_mask; do
mkdir -p $dir/mask/$mask
cp $maskdir/$mask.nii $dir/mask/$mask/
for extract_from_con in 07_any_dist_stay
do echo $extract_from_con 
ix=${extract_from_con:0:2} 
mkdir -p $dir/mask/${mask}/${extract_from_con} 
for subj in {01..25}; 
do for sess in 2; 
do echo $subj $sess 
dataDir=$dir/session_$sess/${extract_from_con} 
fslmeants -i $dataDir/${subj}_wcon_00$ix -m $dir/mask/$mask/$mask -o $dir/mask/$mask/${extract_from_con}/${subj}_${sess}_${mask}_${extract_from_con}.txt 
cp $maskdir/$mask.nii $dir/mask/$mask/${extract_from_con}/ 
done 
done 
done
done

# I think it would be useful to have these ones 03_irrel_dist_stay 02_rel_dist_stay  separately for OneAfterSwitch versus REST
