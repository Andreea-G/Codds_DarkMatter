## Data Analysis for Dark Matter Direct Detection Experiments

### Description 

This program is used to analyze the data produced by Direct Detection (DD) experiments. DD searches look for 
energy deposited within a detector by the collisions between nuclei in a target material and 
Weakly Interactive Massive Particles (WIMPs) belonging to the dark halo of our galaxy. 
The goal of this code is to compare data from various DD experiments and asses the compatibility between 
them in a statistically meaningful way.

The underlying theory and description of the analysis can be found in many papers, a fraction of which are 
listed under the References section below. In particular, this code has been used to produce the Figures 
found in [[1](#1404.7484)-[4](#PhDDisertation)].

At present, it contains the implementation of the analysis for the following types of WIMP interactions:

- Spin-independent (SI) [[1](#1404.7484),[2](#1507.03902)], or spin-dependent (SD) interactions with axial-vector 
(AV) or pseudo-scalar (PS) coupling [[3](#1502.07682)]

- Isospin-conserving or isospin-violating interactions

- Elastic or inelastic scattering (both exothermic and endothermic)

- Contact or long-range interactions (arbitrary mass of the mediator) [[3](#1502.07682)]

- Form factors used: Helm form factor for SI, and the form factors computed by Fitzpatrick et al for SD interactions 
(wherever possible) [[3](#1502.07682),[5,6](#1203.3542)]

Types of analyses implemented:

- Analysis assuming the Standard Halo Model (SHM) [[1](#1404.7484),[2](#1507.03902),[3](#1502.07682)]

- Halo-independent analysis [[1](#1404.7484),[2](#1507.03902),[3](#1502.07682),[7](#1306.5273)]

- Extended maximum likelihood halo-independent analysis (EHI), a recent method for experiments with 
possible dark matter signals and unbinned data [[2](#1507.03902)]


### References
These represent only a few selected references, that contain the description of the analysis implemented in this
code. See also references therein. An initial version of a data analysis code had been implemented in 
Mathematica [[8](#Mathematica_code)].

[1] <a id="1404.7484"></a> G. B. Gelmini, A. Georgescu and J. H. Huh,
  *Direct detection of light Ge-phobic'' exothermic dark matter,*
  JCAP **1407**, 028 (2014) [arXiv:1404.7484 [hep-ph]](http://arxiv.org/abs/1404.7484). 

[2] <a id="1507.03902"></a> G. B. Gelmini, A. Georgescu, P. Gondolo and J. H. Huh,
  *Extended Maximum Likelihood Halo-independent Analysis of Dark Matter Direct Detection Data,*
  [arXiv:1507.03902 [hep-ph]](http://arxiv.org/abs/1507.03902). 

[3] <a id="1502.07682"></a> E. Del Nobile, G. B. Gelmini, A. Georgescu and J. H. Huh,
  *Reevaluation of spin-dependent WIMP-proton interactions as an explanation of the DAMA data,*
  JCAP **1508**, no. 08, 046 (2015) [arXiv:1502.07682 [hep-ph]](http://arxiv.org/abs/1404.7484).
  
[4] <a id="PhDDisertation"></a> A. Georgescu,
 *Tests of WIMP Dark Matter Candidates with Direct Dark Matter Detection Experiments,*
 PhD Dissertation, University of California, Los Angeles (2015) 

[5] <a id="1203.3542"></a> A. L. Fitzpatrick, W. Haxton, E. Katz, N. Lubbers and Y. Xu,
  *The Effective Field Theory of Dark Matter Direct Detection,*
  JCAP **1302** (2013) 004 [arXiv:1203.3542 [hep-ph]](http://arxiv.org/abs/1203.3542). 

[6] A. L. Fitzpatrick, W. Haxton, E. Katz, N. Lubbers and Y. Xu,
DM Direct Detection code written in Mathematica, used to produce the Form Factors described in [[5](#1203.3542)]

[7] <a id="1306.5273"></a> E. Del Nobile, G. Gelmini, P. Gondolo and J. H. Huh,
  *Generalized Halo Independent Comparison of Direct Dark Matter Detection Data,*
  JCAP **1310**, 048 (2013) [arXiv:1306.5273 [hep-ph]](http://arxiv.org/abs/1306.5273). 

[8] <a id="Mathematica_code"</a> A. Georgescu, P. Gondolo, J. H. Huh, 
Data analysis code written in Mathematica [unpublished]


### License

GNU GPL version 3
