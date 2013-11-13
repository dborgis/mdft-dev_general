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
* `read_ck_or_chi` If true the pair correlation functions files specified in 'ck_species' will be read, this is necessary if you want to include any kind of polarization
* `ck_species` specify which type of correlation functions you want to use to describe the solvent, those file are stored in input/direct_correlation_functions... or you can provide your own files in the input directory and use the tag: 'perso'
* `read_chi` If true and if it will read the longitudinal and transverse suceptibily, note that in that case the tag of 'read_ck_or_chi' must be T.
* `polarization` do you want to include any type of polarization: dipolar or multipolar?
* `evaluate_polarization` if dipol will use cdelta and cd to compute the purely dipolar polarization contribution to free energy, if multi will use suceptiblities to compute the multipolar polarozation contribution to free energy
*`temperature` Temperature in Kelvin, it is only involved in KbT terms
*`ref_bulk_density` is the density of the reference fluid, in molecule/Angstrom^3, you should specify the value on the next line.
*`Linearize_entropy`set T if you want to linearize the ideal term of the functunial, if is is the case you need to specify how you want to do that on the next keyword:`if_Linearize_entropy. Usually you do not want to Linearize the ideal term so put the tag F.
*`if_Linearize_entropy` 1: you will linearize wrt n(r) which the density averaged on the angle, if 2: you will linearize wrt rho(r,omega) which is the space and angular dependent density
*`bridge_hard_sphere` T if you want to use an Hard Sphere bridge correction, described in Levesque et al., The Journal of Chemical Physics 137, 034115 (2012), if it is the case, you need also to compute the total free energy for an hard sphere fluid, i.e, the tag on the next keyword `hard_sphere_fluid` MUST ALSO BE T.
*`hard_sphere_fluid` T if you want to treat your solvent as an hard sphere fluid, if it is the case your should also choose an hard sphere as solvent in solvent.in, except if you want to add a Hard sphere correction to a molecular fluid, or if you want to use an HS bridge.(cf supra)
*`hard_sphere_radius`The radius in Angstrom of the Hard sphere solvent is set in line right after this keyword
*`lennard_jones_perturbation_to_hard_spheres`T if you want to add a lennard jones contribution to the hard sphere solvent, this contribution is computed using the WCA theory. If you are not using a HS solvent then the tag MUST BE F. 
### solute.in




### solvent.in
