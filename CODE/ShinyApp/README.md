# A Shiny app for running simple simulations

In this directory, you will find code to run a [Shiny](https://shiny.rstudio.com/) app that allows running simple simulations. The app itself is the file `app.R`, all other files are essentially copies of the ones found in the main code folder.

### Running the app from RStudio

To do this, it suffices to open `app.R` from within RStudio. RStudio detects a Shiny app and therefore shows a "Run App" button, which will run the app.

### Running the app on a Shiny Server

[Shiny Server](https://rstudio.com/products/shiny/shiny-server/) runs on Linux machines and allows to access Shiny apps through a web browser, thereby eliminating the need to run R or RStudio. The content of this folder is sufficient to run the app, so copying this entire folder (and its subfolders) into the Shiny Server root folder is all that is needed, although it could be a good idea to change the name of the folder from `ShinyApp` to something more explicit, since the folder name determines the name of the app on the server.

### Changing neighbourhoods or data file date

We provide a barebones version of the app allowing to run simulations using one of the two neighbourhoods in the paper and with data from 2020-01-28 as in the paper. Using the files and instructions in the main code folder, it is of course possible to change these parameters. 
