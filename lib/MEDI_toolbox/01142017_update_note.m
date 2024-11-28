% This is an updated version of MEDI toolbox.
% 
Important changes:

1. New function to remove readout phase correction - iField_correction().
On some scanner platforms, a readout phase correction is performed that may lead to inconsistencies between echoes within a multiple echo 
acquisition, leading to an unphysical field map when estimated from all echoes together. Similar artifacts are sometimes observed when
bipolar readout is used on certain platforms. iField_correction() removes these phase inconsistencies from the complex data.
Please refer to README.m for details on how to call this function.

Minor changes:

1. Update to Read_Philips_DICOM

Please contact qsmreconstruction@gmail.com with your questions and
suggestion

Regards,
Cornell MRI lab