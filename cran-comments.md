## appac 4.0.3 — initial CRAN submission

This is a new submission; appac is not currently on CRAN. The package
corrects gas-chromatography peak areas for the influence of ambient air
pressure on standard detectors open to the atmosphere (e.g. the
flame-ionization detector).

## Test environments

* Local: Ubuntu, R 4.5.3 (x86_64-pc-linux-gnu)
* GitHub Actions (r-lib/actions/check-r-package):
  * ubuntu-latest, R release and R devel
  * windows-latest, R release
  * macos-latest, R release
* R-hub v2 (rhub::rhub_check), default platform set

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE

  ```
  Maintainer: 'Ruediger Forster <meticulous.measurements@gmail.com>'
  New submission
  ```

  Expected for a first submission.

  The same check flags one URL as possibly invalid:

  ```
  URL: https://www.agilent.com/cs/library/technicaloverviews/public/5989-3425EN.pdf
    From: inst/CITATION
    Status: 403
  ```

  The URL is correct and opens normally in a web browser. The Agilent
  server returns 403 to non-interactive (robot) requests only. It is the
  canonical manufacturer technical overview cited as the method reference
  and is intentionally retained.

## Downstream dependencies

There are currently no downstream dependencies (new package).
