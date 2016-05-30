# MDFT's future

The Molecular Density Functional Theory, by

## Authors

|  |  |  |
| --- | --- | --- |
| Daniel Borgis | PI | Maison de la Simulation
| Maximilien Levesque | PI | École Normale Supérieure
| Guillaume Jeanmairet | PhD student (2011-2014) | now at Ali Alavi's group, Max Planck Institute, Stuttgart
| Ljiljana Stojanovitch | postdoctoral researcher (2013-2014) | 
| Volodymyr Sergiievskyi | postdoctoral researcher (2013-2014) |
| Lu Ding | PhD student (2013-) | works on the resolution of the molecular Ornstein-Zernike equations as projected on a basis of rotational invariants
| Cédric Gageat | PhD student (2015-) | works on core developments and biomolecular applications

## Our git workflow

[A successful git branching model](http://nvie.com/posts/a-successful-git-branching-model/)

## How to Report Bugs Effectively

[A blogpost by Simon Tatham](http://www.chiark.greenend.org.uk/~sgtatham/bugs.html)

## Requirements

- Cmake
- FFTW3

`Cmake` and `FFTW3` can often be install through *apt*, *yum*, *pacman*, *conda*, etc.

## Installation procedure

```sh
git clone https://github.com/maxlevesque/mdft-dev
cd mdft-dev && mkdir build && cd build
cmake ..
make -j
ln -s ../input
```
