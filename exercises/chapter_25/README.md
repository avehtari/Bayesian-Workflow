# Files

* `pga-tour-2024.rds`: Anonymized data for the 2024 PGA Tour golf season. Data courtesy of PGA TOUR ShotLink System.
* `StrokeExportDefinitions.pdf`: Data dictionary for the golf data. Data courtesy of PGA TOUR ShotLink System.

The following R code can be used to construct a smaller dataset that only includes putts and variables to use in the putting model:

```r
pga_tour <- readRDS("pga-tour-2024.rds")
keep_cols <- c(
  "date",
  "player #",
  "course #",
  "hole",
  "distance to pin",
  "1st putt flag",
  "in the hole flag"
)
is_putt <- !is.na(pga_tour[["shot type(s/p/d)"]]) & pga_tour[["shot type(s/p/d)"]] == "P"
putts <- pga_tour[is_putt, keep_cols, drop = FALSE]
```
