#!/bin/sh

## go where this script lives:
cd $(dirname $(readlink -e $BASH_SOURCE))

R -e "rmarkdown::render('celltypeid.Rmd', 'all')"
