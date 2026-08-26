# APRScenario — where things stand

Last updated: 2026-08-26 (submitted to CRAN). Kept in `checks/`, which is in `.Rbuildignore`, so it
never ships to CRAN. This file exists so the state travels with the repo rather
than living on one machine.

## Current position

- `master` = **0.0.4.2**, pushed (latest `41715f7`).
- `tar-package` = `a68f858`, now carries `APRScenario_0.0.4.2.tar.gz`
  alongside the older tarballs.
- **0.0.4.2 was submitted to CRAN on 2026-08-26.** Awaiting the outcome: the
  maintainer confirmation email must be clicked for the submission to go live,
  then CRAN either auto-accepts or a volunteer replies. If they come back with
  questions, the reasoning is all in this file and in `cran-comments.md`.
- CRAN deadline: **2026-09-11**. Miss it and the package is archived.

### CRAN dashboard is churning (still checking published 0.0.4.0)

ERROR count by day: 08-21 two (fedora-clang, fedora-gcc) -> 08-22 three
(+ debian-gcc) -> 08-24 two (r-release-linux, r-oldrel-windows) -> 08-26 one
(r-release-linux only). `bsvarSIGNs` stayed at 3.0 throughout, so this is CRAN
machines churning, not a dependency fix. Do not wait it out: one ERROR on
2026-09-11 is enough to archive the package, and the bug is real for any user
who has bsvarSIGNs 3.0 regardless of the table's colour.

## What was wrong

CRAN reported an ERROR while re-building the vignette:

```
no applicable method for 'forecast' applied to an object of class
"c('PosteriorBSVARSIGN', 'R6')"
```

`bsvarSIGNs` 3.0 (published 2026-08-21) moved which generic owns `forecast`:

| | 2.0 | 3.0 |
|---|---|---|
| NAMESPACE | `importFrom(bsvars, forecast)` | `importFrom(generics, forecast)` |
| method registered in | **bsvars'** S3 table | **generics'** S3 table |
| exports `forecast`? | no | yes |

`bsvars` is still 3.2 and defines its own `forecast` generic, so under 3.0 the
vignette's `bsvars::forecast()` call found no method. A flavour starts failing
as soon as its check machine installs 3.0 — which is why it spread from 2
flavours on 08-21 to 3 on 08-22 while others stayed OK.

This was never a defect in APRScenario's own code: nothing in `R/` or `src/`
calls `forecast()`.

## What was done

Three deletions in `vignettes/APRScenario.Rmd` (mirrored in the gitignored
`doc/` copy), nothing added:

1. The `bsvars::forecast()` call. `f_bvar` was assigned and never read; the
   vignette's unconditional forecast comes from `forc_h()`. Removing it means
   the vignette no longer depends on which package owns the generic, so no
   minimum `bsvarSIGNs` version is needed and it works under both 2.0 and 3.0.
2. `dev.new()` + `grid::grid.newpage()` — opened a graphics device during
   checks, and diverted the IRF panel away from the vignette.
3. `Sys.setlocale("LC_TIME", "English.UTF-8")` — a Windows-only locale name
   that failed silently elsewhere, and mutating locale during a check is
   something CRAN dislikes.

Plus: version bump, NEWS entry, `cran-comments.md` rewritten, `^checks$` added
to `.Rbuildignore`.

## Verification (all re-runnable)

- `Rscript checks/verify_bsvarSIGNs30.R` — installs `bsvarSIGNs` 3.0 into
  `checks/lib30/` (a side library; the system-wide 2.0 install is untouched),
  reproduces the original error, confirms the fix, and re-builds the whole
  vignette under 3.0. Downloads the 3.0 tarball itself if absent.
- `R CMD build . && R CMD check --as-cran --run-donttest -o checks/RCMDCHECK
  APRScenario_0.0.4.2.tar.gz` — last run: **0 ERRORs, 0 WARNINGs, 2 NOTEs**,
  tarball 60,709 bytes, `re-building of vignette outputs ... OK`.

The two NOTEs are local-environment artifacts, both already explained in
`cran-comments.md`: "unable to verify current time", and the non-portable flag
`-mno-omit-leaf-frame-pointer` injected by the local R `Makeconf`.

## What is left

1. ~~Submit to CRAN~~ — **done 2026-08-26**. Watch for the CRAN reply.
2. **Reply to PhoenixEmik** on issue #2 / PR #3 — draft text below under
   Attribution. Needs GitHub write; `gh` is not installed on the Linux box.
3. **Optional, and I would skip it:** `devtools::check_win_devel(".")`. It
   checks against R-devel *the language*, but win-builder is Windows and
   Windows r-devel is green, so it will probably install `bsvarSIGNs` 2.0 and
   never exercise the failing combination. `checks/verify_bsvarSIGNs30.R` is
   the stronger evidence and already passes.

**Done:** `tar-package` was updated on 2026-08-26 (`a68f858`). That branch is
refreshed from master as work lands, not only at release — its own history
("Refresh 0.0.4.0 tarball from master ...") shows the convention.

## Attribution

PhoenixEmik (GitHub) reported the breakage in issue #2 and independently
reached the same diagnosis, then opened PR #3 with a working fix that selected
the generic at runtime. PR #3 was **merged** (`16c1719`) so their commits are in
the history; 0.0.4.2 then removed the call altogether, which made their shim
unnecessary. They are credited in `NEWS.md` under 0.0.4.2 Acknowledgements.

### Draft reply to post on issue #2

> Thanks for the report and the fix — you had the diagnosis right, and fast.
>
> I've merged the PR, then followed up with a commit that removes the
> `forecast()` call altogether: `f_bvar` was assigned but never used, and the
> vignette's unconditional forecast already comes from `forc_h()`. That makes
> the version-dispatch helper unnecessary, so the vignette no longer depends on
> which package owns the generic at all. Released as 0.0.4.2 and submitted to
> CRAN; you're credited in NEWS.md.
>
> For anyone finding this later: the cause is that bsvarSIGNs 3.0 registers
> `forecast.PosteriorBSVARSIGN` on `generics::forecast` rather than on
> `bsvars::forecast`, so a qualified `bsvars::forecast()` call stops
> dispatching. It surfaces on a given CRAN flavour as soon as that machine
> installs 3.0 — which is why the affected flavours shifted from day to day
> rather than staying put.

Note: Tomasz Woźniak (bsvarSIGNs author) advised on issue #2 that this would
resolve itself and to wait. The evidence did not support that — the failure
count rose from 2 to 3 flavours overnight, and the breakage reproduces locally
with `bsvars` 3.2 + `bsvarSIGNs` 3.0.

## Gotchas

- `checks/` must stay in `.Rbuildignore` (`^checks$`). Without it the side
  library ships inside the source tarball: 8 MB instead of 60 KB, plus four
  spurious NOTEs. `checks/apply_0_0_4_2.sh` reverts the working tree, so it
  re-adds this line itself — that was a bug found the hard way.
- `.gitignore` keeps `checks/` artifacts out of git while tracking the scripts.
  Never `git add -A` here without looking: `lib30/` and `RCMDCHECK/` are ~32 MB.
- If you are on a different machine or architecture, delete `checks/lib30`
  before running the verify script — it holds a compiled `.so` built elsewhere,
  and the script only rebuilds when that directory is absent.
- 0.0.4.0 shipped to CRAN, so its `NEWS.md` section must never be relabelled.
  PR #3 did exactly that and it was corrected in 0.0.4.2.
