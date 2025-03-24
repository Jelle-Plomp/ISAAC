# Loading HFR data for analysis
This subdirectory contains some examples to (batch) load (a selection of)
the files containing reconstructed HFR data. The loaded data can then be 
used for further analysis (e.g. Particle Image Velocimetry (PIV)).

Furthermore, code is provided to perform the steps as described in 
[Attenuation Correction and Normalisation for Quantification
of Contrast Enhancement in Ultrasound Images of Carotid Arteries](https://doi.org/10.1016/j.ultrasmedbio.2015.02.010).
by Cheung et al. This can be applied on HFR data to homogenise the image intensity
throughout the image. The implementation is a best matching implementation with
the description in the paper. Implementing this correction by Cheung et al. did not improve
the accuracy of PIV in shadow regions, as shown in the supplementals of our paper.

The PIV analysis itself is not included, since this was performed using
software obtained from Voorneveld et al (see e.g. [High-Frame-Rate Contrast-Enhanced Ultrasound for Velocimetry in the Human Abdominal Aorta](https://doi.org/10.1109/TUFFC.2018.2846416). 


