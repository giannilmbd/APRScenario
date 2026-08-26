# APRScenario 0.0.4.2 — resubmission comments

## Reason for this release

This release fixes the ERROR shown on the package check page, raised while
re-building the vignette:

```
no applicable method for 'forecast' applied to an object of class
"c('PosteriorBSVARSIGN', 'R6')"
```

The cause is a change in the suggested package `bsvarSIGNs`. Up to version 2.0
it registered `forecast.PosteriorBSVARSIGN` on the generic imported from
`bsvars`; version 3.0 registers it on the generic imported from `generics` and
re-exports it. The vignette called `bsvars::forecast()`, which therefore no
longer dispatches. A flavour begins failing as soon as its check machine
installs `bsvarSIGNs` 3.0, which is why the set of affected flavours has
shifted from day to day.

The call has been removed. Its result was never used — the vignette computes
its unconditional forecast with this package's `forc_h()` — and no exported
function of APRScenario calls `forecast()`. The vignette therefore no longer
depends on which package owns the generic, and no minimum version of the
suggested package is imposed. This was verified by installing `bsvarSIGNs` 3.0
into a separate library and re-building the vignette in full: it builds under
both `bsvarSIGNs` 2.0 and 3.0.

Two calls unnecessary to demonstrating the package were also removed from the
vignette: a `dev.new()` that opened a graphics device during checking, and a
`Sys.setlocale("LC_TIME", "English.UTF-8")` using a Windows-only locale name.

## Test environments

* Local: Ubuntu 24.04.4 LTS, R 4.3.3, gcc 13.3.0,
  `R CMD check --as-cran --run-donttest`

## R CMD check results

0 ERRORs, 0 WARNINGs, 2 NOTEs.

Both NOTEs are specific to the local check machine:

* "unable to verify current time" — no network time service is available
  there.
* "Compilation used the following non-portable flag(s):
  -mno-omit-leaf-frame-pointer" — this flag is injected by the local R
  `Makeconf`, not by the package's `src/Makevars`.

## Downstream dependencies

None (no packages depend on APRScenario).
