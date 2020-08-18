# fROCSA
Bruker pulse programs and CSA fitting example for fROCSA 

fROCSA is a solid-state NMR pulse sequence for recoupling chemical shift anisotropy in an indirect dimension during a MAS experiment. Please find more information about this work in our [paper](https://www.biorxiv.org/content/10.1101/2020.07.02.184770v1) on bioaxiv.  A peer-reviewed version of this manuscript is in press at the Journal of Chemical Physics. 

The pulse programs were tested on Bruker Neo spectrometers with Topspin 4.0.7.  We are happy to help you get things working, but use the programs at your own risk.  The anisotropic spectrum must be reversed during processing. 

The needed python packages for the fitting functions can be installed with: `pip install requirements.txt -r `. Please remember to use the correct scaling factors from our paper.

The Mathematica notebooks we wrote to calculate the scaling factors of the fROCSA sequence are also given.

### Citation
*Fritzsching, K. J., *Keeler, E. G., He, C., & McDermott, A. E. (2020). Scaled Recoupling of Chemical Shift Anisotropies at High Magnetic Fields under MAS with Interspersed C-elements. In bioRxiv. https://doi.org/10.1101/2020.07.02.184770
