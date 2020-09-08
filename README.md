# fROCSA
Bruker pulse programs and CSA fitting example for fROCSA 

fROCSA is a solid-state NMR pulse sequence for recoupling chemical shift anisotropy in an indirect dimension during a MAS experiment. Please find more information about this work in our [paper](https://doi.org/10.1063/5.0020682).  A peer-reviewed version of this manuscript is in press at the Journal of Chemical Physics. 

The pulse programs were tested on Bruker Neo spectrometers with Topspin 4.0.7. The anisotropic spectrum must be reversed during processing. 

The needed python packages for the fitting functions can be installed with: `pip install requirements.txt -r `. Please remember to use the correct scaling factors from our paper.

The Mathematica notebooks we wrote to calculate the scaling factors of the fROCSA sequence are also given.

### Citation

*Fritzsching, K. J., *Keeler, E. G., He, C., & McDermott, A. E. (2020). Scaled recoupling of chemical shift anisotropies at high magnetic fields under MAS with interspersed C -elements. The Journal of Chemical Physics, 153(10), 104201. https://doi.org/10.1063/5.0020682

A preprint of this paper is avalable on [BioRxiv](https://doi.org/10.1101/2020.07.02.184770)
