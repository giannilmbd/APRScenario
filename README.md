# APRScenario

This package is a stab implementing [Structural scenario analysis with SVARs](https://doi.org/10.1016/j.jmoneco.2020.04.009) by Juan Antolín-Díaz, Ivan Petrella and Juan F. Rubio-Ramírez, JME (2021) in R

It depends on the output of the package [bsvarSIGNs](https://github.com/bsvars/bsvarSIGNs) although it could be adapted to the output of many other packages by changing the function get_mats.R.

See the application code APRScenario.Rmd. 

## Installation

```
devtools::install_github("giannilmbd/APRScenario", ref = "master")
```
See also the forked bsvarSIGNs with parallelized draws of the rotation matrix
```
devtools::install_github("giannilmbd/bsvarSIGNs", ref = "master")
```

You can also download the tar.gz package from the *tar-package branch* and install in in R with 

```
install.packages("<path>APRScenario_XXXX.tar.gz", repos = NULL, type = "source");
```

For XXXX see the latest version
Use APRScenario.Rmd as a template for the application. 

Feedbacks are welcome.
