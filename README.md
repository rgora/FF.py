ff.py
=====

FF.py is an open source script, written in Python, for preparing a set of input files for finite-field calculations of static electric dipole properties and subsequent estimation of these properties based on the calculated results.
Capabilities

Currently FF.py is capable of preparing the input files for the following quantum chemistry packages:
- Gaussian'09,
- GAMESS(US),
- Molcas.

Unfortunately, the default precision of energy and dipole moments printout in all the above programs must be increased in order to obtain a meaningful results of the first and second hyperpolarizabilities. In the case of Gaussian package the data are extracted from the checkpoint files while for the remaining packages the log files are required.

## Download & Install ##

FF.py is freely available for download from our git repository under the terms of the GNU General Public License. Please note that it requires python (tested on version 2.6.6) and numpy (tested on version 1.4.1).

## Support ##

Please contact the author [mailto:robert.gora@pwr.edu.pl Robert Góra].

## Tutorials & References ##

For the details of finite field approach I can recommend a review by Kurtz and Dudis [1], whereas the details of the iterative Generalized Romberg differentiation procedure can be found in the paper by Medved et al. [2]. This script has been used also in our recent papers [3,4].

The usage of the script itself is rather self-explanatory, however, I have prepared a short tutorial which is available through project's wiki.

## References ##

1. Kurtz, H. A.; Dudis, D. S. In Reviews in Computational Chemistry; Lipkowitz, K. B., Boyd, D. B., Eds.; Wiley: New York, 2007; pp 241-279
2. Medved, M.; Stachová, M.; Jacquemin, D.; André, J.; Perpète, E. A.; J. Mol. Struct. THEOCHEM, 2007, vol. 847, pp 39-46
3. Góra, R. W.; Zaleśny, R.; Zawada, A.; Bartkowiak, W.; Skwara, B.; Papadopoulos, M. G.; Silva, D. L. J. Phys. Chem. A, 2011, 115, 4691-4700