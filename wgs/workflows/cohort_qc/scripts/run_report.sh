#!/bin/bash


DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"



R -e "rmarkdown::render('$DIR/report.Rmd',output_file = '$1',  params=list(label='$2', oncoplot='$3', somatic_plot='$4', summary='$5', burden_plot='$6'))"

