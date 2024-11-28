function Mask = genLiverMask()
rat = load_nii('liver_mask.nii');
Mask = rat.img;
Mask = flip(Mask,1);
Mask = logical(flip(Mask,3));
end