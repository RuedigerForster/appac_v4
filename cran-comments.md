## appac 4.0.4 — maintenance update

This is a documentation-only patch of appac 4.0.3 (currently on CRAN). It
corrects a factual overstatement in the package vignette: the magnitude of the
atmospheric-pressure effect on FID sensitivity was given as "±5 %" (implying a
± band, i.e. roughly twice the real effect); it is up to about 5 % across the
naturally occurring range of ambient pressure. There are no code changes, no
change to the exported API, and no change to dependencies.

## Test environments

* Local: Ubuntu, R 4.5.3 (x86_64-pc-linux-gnu)
* GitHub Actions (r-lib/actions/check-r-package):
  * ubuntu-latest, R release and R devel
  * windows-latest, R release
  * macos-latest, R release
* R-hub v2 (rhub::rhub_check), default platform set

## R CMD check results

0 errors | 0 warnings | 1 note

The check still reports possibly misspelled words in the Description:

```
Bocek (11:15), Clardy (13:11), Janak (11:30), Novak (11:22),
Zeileis (20:56)
```

These are all correctly-spelled author surnames of the cited references
(Boček, Novák and Janák 1969; Ayers and Clardy 1985; Zeileis and others
2002), given in ASCII form. Not misspellings.

The same check flags one URL as possibly invalid:

```
URL: https://www.agilent.com/cs/library/technicaloverviews/public/5989-3425EN.pdf
  From: inst/CITATION
  Status: 403
```

The URL is correct and opens normally in a web browser. The Agilent server
returns 403 to non-interactive (robot) requests only. It is the canonical
manufacturer technical overview cited as the method reference and is
intentionally retained.

## Downstream dependencies

There are currently no downstream dependencies.
