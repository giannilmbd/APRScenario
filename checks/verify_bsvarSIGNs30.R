# Verify APRScenario's vignette against bsvarSIGNs 3.0.
#
# Context: CRAN flagged an ERROR on the two Fedora flavors (2026-08-21):
#   no applicable method for 'forecast' applied to an object of class
#   "c('PosteriorBSVARSIGN', 'R6')"
# bsvarSIGNs 3.0 moved its forecast() method from the bsvars generic to the
# generics generic. Those builders had 3.0; every other flavor (and this
# machine) still had 2.0, which is why local checks passed.
#
# The Fedora check aborted at the forecast call, so nothing downstream of it
# has ever run under 3.0. This script installs 3.0 into a SIDE LIBRARY (the
# system-wide 2.0 install is left untouched) and renders the full vignette,
# to confirm the rest of the vignette works under 3.0 too.
#
# Run:  Rscript checks/verify_bsvarSIGNs30.R
# Log:  checks/verify_bsvarSIGNs30.log

# run from the project root
proj <- getwd()
if (!dir.exists(file.path(proj, "vignettes")))
  stop("run this from the APRScenario project root")

lib30   <- file.path(proj, "checks", "lib30")
tarball <- file.path(proj, "checks", "bsvarSIGNs_3.0.tar.gz")
vig     <- file.path(proj, "vignettes", "APRScenario.Rmd")

dir.create(lib30, showWarnings = FALSE, recursive = TRUE)

# --- 1. install bsvarSIGNs 3.0 into the side library -----------------------
if (!dir.exists(file.path(lib30, "bsvarSIGNs"))) {
  message("Installing bsvarSIGNs 3.0 into ", lib30, " ...")
  install.packages(tarball, lib = lib30, repos = NULL, type = "source",
                   INSTALL_opts = "--no-multiarch")
}

# side library FIRST so 3.0 shadows the system-wide 2.0
.libPaths(c(lib30, .libPaths()))

stopifnot(
  "bsvarSIGNs 3.0 not picked up" =
    as.character(packageVersion("bsvarSIGNs")) == "3.0"
)
message("bsvarSIGNs: ", packageVersion("bsvarSIGNs"),
        "   bsvars: ",  packageVersion("bsvars"),
        "   APRScenario: ", packageVersion("APRScenario"))

# --- 2. confirm the reported CRAN failure, and that it is generic-ownership -
data(NKdata, package = "APRScenario")
spec <- bsvarSIGNs::specify_bsvarSIGN$new(
  data = as.matrix(NKdata[, 2:4]), p = 1)
post <- bsvars::estimate(spec, S = 10)

message("class(post): ", paste(class(post), collapse = ", "))

old_generic <- tryCatch({ bsvars::forecast(post, horizon = 3); "dispatched" },
                        error = function(e) conditionMessage(e))
new_generic <- tryCatch({ bsvarSIGNs::forecast(post, horizon = 3); "dispatched" },
                        error = function(e) conditionMessage(e))
message("bsvars::forecast()     -> ", old_generic)
message("bsvarSIGNs::forecast() -> ", new_generic)

# --- 3. confirm the API APRScenario actually consumes is unchanged in 3.0 ---
stopifnot(
  "posterior$posterior$A missing"      = !is.null(post$posterior$A),
  "posterior$posterior$B missing"      = !is.null(post$posterior$B),
  "specification$data_matrices missing"= !is.null(spec$data_matrices$Y)
)
mats <- APRScenario::gen_mats(posterior = post, specification = spec)
message("gen_mats() OK: n_var=", mats$n_var, " n_p=", mats$n_p,
        " dim(M)=", paste(dim(mats$M), collapse = "x"))

# --- 4. render the whole vignette under 3.0 --------------------------------
outdir <- file.path(proj, "checks", "vignette_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
message("Rendering vignette under bsvarSIGNs 3.0 ...")
rmarkdown::render(vig, output_dir = outdir, quiet = FALSE,
                  envir = new.env(parent = globalenv()))
message("VIGNETTE BUILD OK under bsvarSIGNs 3.0")
