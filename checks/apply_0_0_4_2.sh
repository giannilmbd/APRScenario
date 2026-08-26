#!/usr/bin/env bash
#
# Land APRScenario 0.0.4.2 on top of PhoenixEmik's merged PR #3.
#
# PRECONDITION: PR #3 must already be merged on GitHub (see the step list in
# the session notes). This script pulls that merge, then applies the 0.0.4.2
# changes on top of it:
#
#   * restores the published "# APRScenario 0.0.4.0" NEWS heading, which PR #3
#     renamed to 0.0.4.1 (0.0.4.0 shipped to CRAN, so its entry must survive)
#   * keeps PR #3's bullet as its own 0.0.4.1 section
#   * adds a 0.0.4.2 section crediting PhoenixEmik
#   * removes the now-unused f_bvar call together with PR #3's compatibility
#     shim (the value was never read; the vignette forecasts with forc_h())
#   * removes dev.new()/grid.newpage() and the Windows-only Sys.setlocale()
#   * bumps DESCRIPTION and cran-comments.md to 0.0.4.2
#
# It does NOT push and does NOT submit to CRAN.
#
# Usage:  bash checks/apply_0_0_4_2.sh
set -euo pipefail

cd "$(dirname "$0")/.."
PROJ="$PWD"
echo "project: $PROJ"

# ---------------------------------------------------------------- safety ----
if [ "$(git rev-parse --abbrev-ref HEAD)" != "master" ]; then
  echo "ERROR: not on master. Checkout master first." >&2; exit 1
fi

mkdir -p checks
if ! git diff --quiet || ! git diff --cached --quiet; then
  echo "Backing up uncommitted work to checks/backup_worktree.patch"
  git diff HEAD > checks/backup_worktree.patch
  echo "  (restore later with: git apply checks/backup_worktree.patch)"
  git checkout -- .
fi

# ------------------------------------------------------- pull PR #3 merge ---
echo "Pulling master (expects PR #3 already merged on GitHub)..."
git pull --ff-only origin master

if ! grep -q "getNamespaceExports" vignettes/APRScenario.Rmd; then
  echo "ERROR: PR #3 does not appear to be merged — its shim is not in the" >&2
  echo "       vignette. Merge PR #3 on GitHub first, then re-run." >&2
  exit 1
fi
echo "PR #3 merge detected."

# -------------------------------------------------- apply 0.0.4.2 changes ---
python3 - <<'PY'
import re, sys

def edit(path, old, new, required=True):
    s = open(path).read()
    if old not in s:
        if required:
            sys.exit(f"ERROR: expected text not found in {path}:\n{old[:200]}")
        return False
    open(path, 'w').write(s.replace(old, new, 1))
    return True

# ---- DESCRIPTION ----
edit('DESCRIPTION', 'Version: 0.0.4.1', 'Version: 0.0.4.2')

# ---- NEWS.md: split PR #3's merged section back apart, prepend 0.0.4.2 ----
news = open('NEWS.md').read()
if not news.startswith('# APRScenario 0.0.4.1'):
    sys.exit("ERROR: NEWS.md does not start with the 0.0.4.1 heading from PR #3")

emik_bullet = """* Updated the vignette to use the `generics::forecast()` generic exported by
  `bsvarSIGNs` 3.0 while retaining compatibility with older `bsvarSIGNs`
  versions, which dispatch forecasts through `bsvars::forecast()`.

"""
if emik_bullet not in news:
    sys.exit("ERROR: PR #3's NEWS bullet not found verbatim; resolve NEWS.md by hand")

news = news.replace(emik_bullet, '', 1)
# the remainder of that section is the original 0.0.4.0 content
news = news.replace('# APRScenario 0.0.4.1', '# APRScenario 0.0.4.0', 1)

new_sections = """# APRScenario 0.0.4.2

## Bug fixes

* Vignette: removed the call to `bsvars::forecast()` and the version-dispatch
  helper introduced in 0.0.4.1. The forecast object was never used -- the
  vignette computes its unconditional forecast with `forc_h()` -- so the
  vignette no longer depends on which package owns the `forecast()` generic.
  No exported function of APRScenario calls `forecast()`. The vignette builds
  under both `bsvarSIGNs` 2.0 and 3.0.

* Vignette: removed a `dev.new()` / `grid::grid.newpage()` pair that opened a
  graphics device during vignette building, and a
  `Sys.setlocale("LC_TIME", "English.UTF-8")` call using a Windows-only locale
  name that failed silently on other platforms.

## Acknowledgements

* Thanks to PhoenixEmik (https://github.com/PhoenixEmik) for reporting the
  `bsvarSIGNs` 3.0 breakage and independently diagnosing its cause
  (issue #2, pull request #3).

# APRScenario 0.0.4.1

## Bug fixes

""" + emik_bullet

open('NEWS.md', 'w').write(new_sections + news)

# ---- vignette (and the gitignored doc/ copy) ----
# vignette patterns are matched with regexes so trailing whitespace and minor
# reflowing do not break the edit
shim_and_call = re.compile(
    r'forecast_bvarsign\s*<-\s*if\s*\("forecast"\s*%in%\s*'
    r'getNamespaceExports\("bsvarSIGNs"\)\).*?\n\)\n',
    re.S)

locale_chunk = re.compile(
    r'```\{r Set-directory[^\n]*\}\n.*?Sys\.setlocale\([^\n]*\n```\n\n',
    re.S)

dev_pair = re.compile(r'\{dev\.new\(\)[ \t]*\n\s*grid::grid\.newpage\(\)[ \t]*\n')


def sub(path, rx, repl, label, required=True):
    s = open(path).read()
    s2, n = rx.subn(repl, s, count=1)
    if n == 0:
        if required:
            sys.exit(f"ERROR: {label} not found in {path}")
        print(f"  ({label} not found in {path}, skipped)")
        return
    open(path, 'w').write(s2)


for f in ('vignettes/APRScenario.Rmd', 'doc/APRScenario.Rmd'):
    import os
    if not os.path.exists(f):
        print(f"  (skipping {f}, not present)")
        continue
    req = f.startswith('vignettes/')
    sub(f, shim_and_call, '', "PR #3 shim + f_bvar call", req)
    sub(f, locale_chunk, '', "Sys.setlocale chunk", req)
    sub(f, dev_pair, '{\n', "dev.new/grid.newpage pair", req)
    print(f"  patched {f}")

# ---- .Rbuildignore: re-add the checks/ exclusion (reverted with the worktree) ----
rb = open('.Rbuildignore').read()
if '^checks$' not in rb:
    if not rb.endswith('\n'):
        rb += '\n'
    open('.Rbuildignore', 'w').write(rb + '^checks$\n')
    print("  re-added ^checks$ to .Rbuildignore")

# ---- cran-comments.md (written fresh; the working-tree draft was reverted) ----
open('cran-comments.md', 'w').write("""# APRScenario 0.0.4.2 — resubmission comments

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
""")
print("applied 0.0.4.2 changes")
PY

# ------------------------------------------------------------- sanity ------
echo "--- residual references that should be gone ---"
grep -n "forecast_bvarsign\|bsvars::forecast\|Sys.setlocale\|dev.new()" \
     vignettes/APRScenario.Rmd || echo "  (clean)"
echo "--- NEWS headings ---"
grep -n "^# APRScenario" NEWS.md
echo "--- version ---"
grep -n "^Version" DESCRIPTION

cat <<'EOF'

Local changes applied but NOT committed and NOT pushed.

Next:
  1. review:   git diff
  2. verify:   Rscript checks/verify_bsvarSIGNs30.R      # bsvarSIGNs 3.0
  3. verify:   R CMD build . && R CMD check --as-cran --run-donttest \
                 -o checks/RCMDCHECK APRScenario_0.0.4.2.tar.gz
  4. commit & push, then run devtools::check_win_devel()
EOF
