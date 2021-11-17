#!/bin/bash -x
## data is now under gerrit:/data/groups/gen/philip/CHETAH/Seurat_objects
## and locally in /Users/philiplijnzaad/data/elixir-scrnaseq-course2019


## symlink as follows:
## ovarian.rds -> /Users/philip/elixir-scrnaseq-course2019/data/ovarian1200.rds
## chetah.ref.rds -> /Users/philip/elixir-scrnaseq-course2019/data/chetah.ref.rds
## ribosomalgenes.rds -> /Users/philip/elixir-scrnaseq-course2019/data/ribosomalgenes.rds
## cellcyclemarkers.rds -> /Users/philip/elixir-scrnaseq-course2019/data/cellcyclemarkers.rds
## chetah.ref.singler.rds -> /Users/philip/elixir-scrnaseq-course2019/data/chetah.ref.singler.rds

## go where this script lives:
cd $(dirname $(readlink -e $BASH_SOURCE))

cd ./session-celltypeid_files/ || exit 55

srcdir=$HOME/data/elixir-scrnaseq-course2019/

for i in $srcdir/*.rds; do
  ln -s $i ./
done
