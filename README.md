# MDFT's future

The Molecular Density Functional Theory, by

- Daniel Borgis, PI, Maison de la Simulation  
- Maximilien Levesque, PI, École Normale Supérieure
- Guillaume Jeanmairet, as a PhD student (2011-2014), worked on the nP functional for water.
- Ljiljana Stojanovitch, as a postdoctoral researcher (2013-2014), worked on reference MD, PMFs and hydrophobicity.
- Volodymyr Sergiievskyi, as a postdoctoral researcher (2013-2014), worked on pressure corrections and hard sphere bridges.
- Lu Ding, as a PhD student (2013-), works on the resolution of the molecular Ornstein-Zernike equations as projected on a basis of rotational invariants.
- Cédric Gageat, as a PhD student (2015-), works on core developments and biomolecular applications.
- Yacine Ould-Rouis, as an engineer in high performance computing (2017-), works on openmp implementation and vectorization.

## Git, github, issues etc.

[A successful git branching model](http://nvie.com/posts/a-successful-git-branching-model/)

[How to reports bugs effectively](http://www.chiark.greenend.org.uk/~sgtatham/bugs.html)

[How to write a git commit message](http://chris.beams.io/posts/git-commit/)


## Requirements

- Cmake

`Cmake` can often be install through *apt*, *yum*, *pacman*, *conda*, etc.  with `sudo apt install cmake`.

- gfortran

We also recommend version 5 or 6 of gfortran, even if some late versions 4.x should work.  
To install gfortran 6 in Ubuntu (16.04), you need to add the dedicated repository and use apt:  
```sh
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt update
sudo apt install gfortran-6
```

## Installation procedure

To proceed with the compilation of mdft-dev:

```sh
git clone --recursive https://github.com/maxlevesque/mdft-dev
cd mdft-dev
mkdir build
cd build
cmake ..
make -j 4
```

If you want a developer version, with all warnings and debugging flags on, change `cmake ..` to `cmake -DCMAKE_BUILD_TYPE=DEBUG ..`.  There are lots of warnings during the compilation of mdft-dev right now.

mdft-dev requires the FFTW3 library. mdft-dev checks for it on default folders on your computer. If it does not find it, it
downloads it from the FFTW3 website and compiles it, which takes 4 minutes on my laptop. Note it is in this last case very verbose since we compile FFTW3 twice.

## Input files

You should still be in the build dir, that is `mdft-dev/build`. There you need lots of input files that should be transparent to you as a user. Two files should nevertheless be created by you: `dft.in` and `solute.in`. Minimal examples are found in  `data/examples/`.  

## Use mdft

From the build dir, a simple `./mdft-dev` will produce some outputs to the terminal (to stdout), create a folder for all output files called `output`.
