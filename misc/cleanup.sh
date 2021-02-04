#!/bin/bash

git_dir=$(git rev-parse --show-toplevel)

rm $git_dir/Thesis-*.cpt
rm $git_dir/Thesis.aux
rm $git_dir/Thesis.bbl
rm $git_dir/Thesis.blg
rm $git_dir/Thesis.fdb_latexmk
rm $git_dir/Thesis.fls
rm $git_dir/Thesis.log
rm $git_dir/Thesis.out
rm $git_dir/Thesis.pdf
rm $git_dir/Thesis.synctex.gz
rm $git_dir/Thesis.toc
rm $git_dir/figures/*eps-converted-to.pdf
