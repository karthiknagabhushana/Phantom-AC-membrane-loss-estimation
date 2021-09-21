# Phantom-AC-membrane-loss-estimation

This code is based on the following papers -  
1) K. Nagabhushana, W. D. O'Brien and A. Han, "Technique to compensate for unknown laminate transmission loss in phantom attenuation measurements," 2020 IEEE International Ultrasonics Symposium (IUS), 2020, pp. 1-4, doi: 10.1109/IUS46767.2020.9251582.
2) K. Nagabhushana, Q. Wang, W. D. O'Brien and A. Han, , "Pulse-echo Technique to Compensate Laminate Membrane Transmission Loss in Phantom Attenuation Measurements" (under preparation)

The papers aims to simulataneously measure the reference phantom AC and membrane loss using a pulse-echo technique. 

The code includes two components - 1) k-Wave simulation of a use-case (phantom, membrane and transducer), 2) Code to process the RF-data from simulations and estimate the phantom material AC and membrane transmission coefficient. This code can also be repurposed to process the experimental RF-data. 

k-Wave MATLAB codes (for 10MHz, F4 transducer) to generate the input HDF5 files for CPU based k-wave execution. Modify the transducer, phantom and membrane properties to simulate a different use case. 
1. Txdc_10M_F4_3D_ref.m - simulates Setup1 from Paper1.
2. Txdc_10M_F4_3D_att.m - simulates Setup2 from Paper1.
3. Txdc_10M_F4_3D_att_top.m - simulates Setup3 from Paper1.

Batch scripts for CPU based k-wave execution
submit_ref.sh, submit_att.sh, submit_att_top.sh

AC_process.m - MATLAB code for processing the RF-data and estimating phantom thickness, phantom material sound speed, membrane transmission loss, and phantom AC. 
  <functions> 
  1. snip_array.m - snips the relevant time-domain echoes from the RF-data. 
  2. spect_atten.m, zfft.m - compute echo spectrum.
  3. speedinwater.m - calculates the water sound speed based on water temperature. This can be used while processing experimental data. 

  

