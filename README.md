# MDFT

## What for

MDFT is classical, molecular density functional theory code. It is aimed at computing density profiles of solvent around a solute
and the free energy of solvation.

## Authors

MDFT is developped in the group of Daniel Borgis at the Ecole Normale Superieure, Paris, France, by
	- Daniel Borgis
	- Maximilien Levesque
	- Guillaume Jeanmairet

## Installation

For now, a simple `make` will do all the necessary stuff.

## Input File

Three files are necessary.
-   `dft.in` contains simulation informations
-   `solute.in` contains informations about the solute (force field, position, etc)
-   `solvent.in` contains all information about the solvent (force field, position, etc)

### dft.in

#### Preliminary information

`dft.in` is read by the program in any order. Lines may be thus be interchanged.
Only usefull tags will be omitted, but no default choice is made for the user.

`MDFT` is looking for tags. For instance, the whole `dft.in` is parsed in order to find `nfft1`. The value of which is then read
and feeds the code.

Comments should be declared by symbol `#`. Everything that follows this symbol is not considered by `MDFT`.
For instance, in `nfft1 = 32 # 64 in the previous calculation`, only `nfft1 = 32` is considered by `MDFT`.

The rule is `tag = info`. There should be a space before and another after the equal sign.

Empty and blank lines are not considered.

* `Lx` Length of the orthorombic supercell in x direction
* `Ly` Length of the orthorombic supercell in y direction
* `Lz` Length of the orthorombic supercell in z direction
* `nfft1` Number of grid points in x direction
* `nfft2` Number of grid points in y direction
* `nfft3` Number of grid points in z direction
* `quadrature` : quadrature to be used for angular integration, i.e., for the discretization of the angular grid.
    - `L` for [Lebedev quadrature](http://en.wikipedia.org/wiki/Lebedev_quadrature) (recommanded)
    - `GL` for [Gauss-Legendre quadrature](http://en.wikipedia.org/wiki/Gaussian_quadrature)
* `order_of_quadrature` Order of the angular integration. Indirectly define the number of discret angles of the angular grid.
    - `3` is recommanded for water
    - For Gauss-Legendre quadratures, the number of angles is 2*2^order.
    - For Lebedev quadratures, it is 2/3 of the number of angles one would have for Gauss-Legendre quadratures. Works only for 6, 14, 26 and 38.
* `nb_psi` Number of discrete angles for third Euler angle, i.e., for the rotation around the molecular axis.
    - `1` for stockmayer and other linear molecules
    - `4` or more are recommanded for water
* `sym_order` Order of symetry of the main axis of the solvent molecule
    - `2` for a water molecule that is of molecular point group C2v
* `nb_implicit_species` Number of solvent species
    - `1` if the solvent is pure
* `mole_fractions` Mole fraction of each solvent species. Sum of mole fractions of all solvent species is 1.
    - `1.0` if the solvent is pure
* `verbose`
    - `T` for lots of information printed to terminal. Advised, e.g., for debugging purpose
    - `F` for production runs
* `read_ck_or_chi` If true the pair correlation functions files specified in 'ck_species' will be read, this is necessary if you want to include any kind of polarization
* `ck_species` specify which type of correlation functions you want to use to describe the solvent, those file are stored in input/direct_correlation_functions... or you can provide your own files in the input directory and use the tag: 'perso'


### solute.in




### solvent.in
