# APRScenario 0.0.4.2 — resubmission comments

## Reason for this release

This release corrects the ERROR shown on the CRAN check page for the R-devel
Linux flavours (`debian-gcc`, `fedora-clang`, `fedora-gcc`), raised while
re-building the vignette:

```
no applicable method for 'forecast' applied to an object of class
"c('PosteriorBSVARSIGN', 'R6')"
```

The cause is a change in the suggested package `bsvarSIGNs`. Up to version 2.0
it registered its `forecast()` method for `PosteriorBSVARSIGN` objects on the
generic imported from `bsvars`; version 3.0 registers the method on the generic
imported from `generics` and re-exports it. The vignette called
`bsvars::forecast()`, which therefore no longer dispatches. The failure appears
on a flavour as soon as that check machine installs `bsvarSIGNs` 3.0.

The call has been removed. Its result was never used — the vignette computes
its unconditional forecast with this package's `forc_h()` — and no exported
function of APRScenario calls `forecast()`, so no user-visible functionality
changes. The vignette builds under both `bsvarSIGNs` 2.0 and 3.0, verified by
installing 3.0 into a separate library and re-building it in full, so no
minimum version of the suggested package is imposed.

Two unrelated calls that are unnecessary during vignette building were also
removed: a `dev.new()` that opened a graphics device, and a
`Sys.setlocale("LC_TIME", "English.UTF-8")` using a Windows-only locale name.

## Test environments

* Local: Ubuntu 24.04.4 LTS, R 4.3.3, gcc 13.3.0,
  `R CMD check --as-cran --run-donttest`

## R CMD check results

0 ERRORs, 0 WARNINGs, 2 NOTEs.

Both NOTEs are specific to the local check machine:

* "unable to verify current time" — no network time service is available on
  the check machine.
* "Compilation used the following non-portable flag(s):
  -mno-omit-leaf-frame-pointer" — this flag is injected by the local R
  `Makeconf`, not by the package's `src/Makevars`.

## Downstream dependencies

None (no packages depend on APRScenario).
