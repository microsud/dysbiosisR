

<!-- badges: start -->  
[![R-CMD-check](https://github.com/microsud/dysbiosisR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/microsud/dysbiosisR/actions/workflows/R-CMD-check.yaml) [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://github.com/microsud/dysbiosisR/blob/master/LICENSE.md) ![GitHub Repo stars](https://img.shields.io/github/stars/microsud/dysbiosisR?style=social) [![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://lifecycle.r-lib.org/articles/stages.html#experimental)   
<!-- badges: end -->

# dysbiosisR <img src="man/figures/logo.png" align="right" height="130"/>

The *dysbiosisR* package is a collection of functions to calculate some common 
dysbiosis measures from microbiota profiling data. The package currently supports 
microbiome data organized as phyloseq objects.  


## Installation

You can install the development version of `dysbiosisR` like this:
```
    devtools::install_github("microsud/dysbiosisR")
``` 

## Load package  
``` 
    library(dysbiosisR)
``` 
The package is currently in the Beta version and a stable version will be made 
available soon. We welcome feedback and suggestions.

## How to use   

Check out the 
[dysbiosisR vignette](https://microsud.github.io/dysbiosisR/articles/Introduction.html) and 
[documentation of functions](https://microsud.github.io/dysbiosisR/reference/index.html) for more 
details.    

## Contribute

Contributions and feedback are very welcome:  

-   [Raise Issue](https://github.com/microsud/dysbiosisR/issues)
-   [Star us on the Github page](https://github.com/microsud/dysbiosisR)

## Cite

Shetty SA, et al., 2022. *dysbiosisR*: An R package to calculate common dysbiosis measures.

### References

1.  Santiago M, et al., 2019. Microbiome predictors of dysbiosis and VRE decolonization in patients with recurrent C. difficile infections in a multi-center retrospective study. AIMS Microbiol 5:1–18. 10.3934/microbiol.2019.1.1.

2.  Lloyd-Price J, Arze C, Ananthakrishnan AN, et al., 2019. Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases. Nature 569:655–662. 10.1038/s41586-019-1237-9.

3.  AlShawaqfeh MK, et al., 2017. A dysbiosis index to assess microbial changes in fecal samples of dogs with chronic inflammatory enteropathy. FEMS Microbiol Ecol 93:1–8. 10.1093/femsec/fix136.

4.  Montassier E, Al-Ghalith GA, Hillmann B, Viskocil K, Kabage AJ, McKinlay CE, Sadowsky MJ, Khoruts A, Knights D. 2018. CLOUD: a non-parametric detection test for microbiome outliers. Microbiome 6:137. 10.1186/s40168-018-0514-4.

5.  Saffouri GB, Shields-Cutler RR, Chen J, Yang Y, Lekatz HR, Hale VL, Cho JM, Battaglioli EJ, Bhattarai Y, Thompson KJ, Kalari KK, Behera G, Berry JC, Peters SA, Patel R, Schuetz AN, Faith JJ, Camilleri M, Sonnenburg JL, Farrugia G, Swann JR, Grover M, Knights D, Kashyap PC. 2019. Small intestinal microbial dysbiosis underlies symptoms associated with functional gastrointestinal disorders. Nat Commun 10:2012. 10.1038/s41467-019-09964-7.

6.  Wei, S., Bahl, M.I., Baunwall, S.M.D., Hvas, C.L. and Licht, T.R., 2021. Determining gut microbial dysbiosis: a review of applied indexes for assessment of intestinal microbiota imbalances. Applied and Environmental Microbiology, 87(11), pp.e00395-21.

## Acknowledgement

This work was done at the Department of Medical Microbiology and Infection Prevention at UMC Groningen in collaboration with the Microbiome Group at the RIVM.

## Contact

Sudarshan A. Shetty (sudarshanshetty9[\@]gmail[dot]com)  
Susana Fuentes (susana[dot]fuentes[\@]rivm[dot]nl)
