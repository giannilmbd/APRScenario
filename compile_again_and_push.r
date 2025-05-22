# recompile everything and upload to github

system('rm src/*.so')
Rcpp::compileAttributes(); devtools::document();devtools::build()
system('git add .;git commit -m "recompiled everything"; git push origin master')
system('cd ../only_tar;mv ../APRScenario_0.0.3.0.tar.gz .; git add -f APRScenario_0.0.3.0.tar.gz; git commit -m "added new tarball"; git push origin tar-package')
print("done")
