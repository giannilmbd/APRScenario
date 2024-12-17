# APRScenario
This package is a stab implementing [Structural scenario analysis with SVARs](https://www.sciencedirect.com/science/article/abs/pii/S0304393220300751) by Juan Antolín-Díaz, Ivan Petrella and Juan F. Rubio-Ramírez, JME (2021) in R

It depends on the output of the package [bsvarSIGNs](https://github.com/bsvars/bsvarSIGNs) although it could be adapted to the output of many other packages by changing the function get_mats.R.

See the application in the quarto (R) code APRScenario.qmd. 

Download the tar.gz file and install in in R with 
```
install.packages("<path>APRScenario_0.0.0.9000.tar.gz", repos = NULL, type = "source");
```
Use APRScenario.qmd as a template for the application. 

Feedbacks are welcome.
