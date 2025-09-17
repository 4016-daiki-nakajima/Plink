# Plink: A Fast Framework for Modeling Sound from Tiny Rigid Bodies

## Installation

Building is handled by CMake:

    git submodule update --init --recursive
    mkdir build && cd build
	cmake -DCMAKE_BUILD_TYPE=Release ..
    make -j4

## Run 

To run an example:

    ./Plink sphere
    
## Dataset

You can download mesh obj files such as obj for bunny and dragon from the [KleinPAT dataset](https://graphics.stanford.edu/projects/kleinpat/kleinpat-dataset/dataset_table.html). 

## Acknowledgement

This project is build on top of the [Simple Modal](https://github.com/zhehaoli1999/SimpleModal) Project by [Zhehao Li](https://zhehaoli1999.github.io/).
