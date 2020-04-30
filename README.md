# macro-design-app
Shiny app for designing large-scale camera trap surveys using balanced acceptance sampling (BAS). The app allows you to upload a shapefile demarcating a (large) study area, and generates spatially balanced points within this survey area. Each of these points is intended to be the rough centroid of a future SECR surveys (e.g. camera trap surveys). Extensions allow for stratification of the study area, and the exclusions of parts that have already been surveyed. The BAS points can be exported as a csv file.

Click this link to run the app in your web browser (there is a limit to the number of concurrent users):

[macro-design-app](https://iandurbach.shinyapps.io/macro-design-app/)

To run locally on your own machine, paste and run this code at the R command prompt:

```r
library(shiny)
runGitHub("macro-design-app", "iandurbach")
```

This requires that packages `shiny`, `dplyr`, `sf`, and `leaflet` are already installed. Otherwise, download or fork the repo as usual. Some test shapefiles are provided in the *data* folder of this repository.

Note the binder link below fails, don't use it (to do).

<!-- badges: start -->
[![Launch Rstudio Binder](http://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/iandurbach/macro-design-app/master?urlpath=rstudio)
<!-- badges: end -->
