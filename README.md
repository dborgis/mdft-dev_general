# MDFT's future

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
