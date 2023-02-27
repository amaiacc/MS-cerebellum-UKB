#-------------------------------------------------------------------------------
# Globals
#-------------------------------------------------------------------------------

rm(list=ls())

# devtools::install_github("rstudio/rmarkdown")

library(knitr)
library(rmarkdown)
# library(pander)
#-------------------------------------
args <- commandArgs(TRUE)
print(args)

#------------------------------------
script_dir="/general_scripts/genotyping/ldsc/"

## general
# script="ldsc_h2_output_UKB_BIG40.R"
# script="ldsc_rg_atlas_output_UKB_BIG40.R"

## cerebellum
script="cerebellum/ldsc_UKB_BIG40_cerebellum_output.R"
# script="cerebellum/ldsc_UKB_BIG40_cerebellum_behav_output.R"
# script="cerebellum/ldsc_UKB_BIG40_cerebellum_cortical_output.R"
# script="cerebellum/ldsc_UKB_BIG40_cerebellum_subcortical_output.R"

#------------------------------------
# Set directories, dependent on  system
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
#------------------------------------
# define file to run
working_dir=paste(dir,script_dir,sep="")
file=paste(working_dir,script,sep="")
#-------------------------------------------------------------------------------
# Main function
#-------------------------------------------------------------------------------

# rmarkdown::render function uses knitr::spin under the hood

opts_chunk$set(include=TRUE, echo=FALSE, error = FALSE, warning = FALSE, message=FALSE, results = "asis", tidy=TRUE) # , fig.width=8, fig.height=6
render(file,
       output_format = "html_document", #"pdf_document",
       output_dir = paste(working_dir,"reports",sep="")
)

