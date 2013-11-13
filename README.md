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
*`hs_functional` which type of functional you want to use for the hard sphere fluid, either Carnahan-Starling(CS) or Percus-Yevick (PY)

Hereafter you described what kind of potential you want to use to compute the external potential.
*`poisson_solver`if T the electrostatic potential generated by the solute will be computed using a poisson solver, i.e using the charge density of the solute, it is quicker than the folowwing way, so it should be used in most of the case.
*`point_charge_electrostatic`if T the electrostatic potential generated by the solute will be computed using a direct way i.e V=q/(4pi.epsilon0.r) this is slower than the previous way.
WARNING: Note that if you set T to the two previous keywords, you will computed the electrostatic potential due to the solute TWICE, this should never be the case.

*`reuse_density`if T instead of starting the minimization from a guess density, it will use as a starting point the density stored in density.bin, so if you want to restart a minimization from a previous one you can move the file output/density.in into input/ and set T to this keyword
*`hard_sphere_solute` if T it will add an Hard Sphere in the center of the box, with a radius that you can set in the next line. If you want to use a purely hard sphere for solute , you also need to use the solute called "Hard Sphere" in solute.in

*`HARD WALL` the number of hard wall you want to introduce in the supercell, those walls are described by the n next lines (n beeing the number of hard walls you decided to introduce). There are described by 7 reals, the fist is the thickness of the wall, the next 3 are the coordinates of a vector normal to the surface of the wall, the next 3 are coordinates on Ox Oy Oz axis of a point belonging the median plane of the wall.
*`vdw_hard_core` ???
*`hard_cylinder` if T it will include an hard cylinder at the center of the box the axis of the cylinder in on the z direction
*`radius_of_hard_cylinder`here you specify what radius you want for the hard cylinder, in Angstrom
*`personal_vext` if you want to use a specific form for the external potential that is not possible to use with the descriptors included in this dft.in file, you can set T to this keyword and then edit the file src/personal_vext.f90 where you can compute your specific external potential
*`vext_hard_square_well` T if you want to use a square well potential as external potential
*`hard_well_total_length`specify the length of the square well.  potential is zero from Lz/2-hard_well_total_length/2 to Lz/2+hard_well_total_length/2 and \infty elsewhere
*`purely_repulsive_solute`if T it will add the repulsive part of a LJ to your solute, as described in Dzubiella et al., J. Chem. Phys., Vol. 121, No. 11, 15 September 2004. If you want to include only this repulsive part, you need to use the solute called "Hard Sphere" in solute.in
*`radius_of_purely_repulsive_solute` is the radius of the purely repulsive sphere as descibed in Dzubiella et al., J. Chem. Phys., Vol. 121, No. 11, 15 September 2004.

*`minimizer` The minimizer you want to use, bfgs is the more convincing, this should be you 1st choice
*`maximum_iteration_nbr` is the number of iterations afterwich you will stop the minimization if a minimum has not been found yet. Usually if you have not converged after 50 iterations you can give up, something went wrong.
*`epsg` is one criterium of convergence, if at iteration n (F(n-1)-F(n))/F(n)< this value, the minimizer will consider that a minimum is reach, by default we use 0.0001
*`pgtol`is the other criterium of convergence if the sum of the norm of dF[rho(r,Omega)]/drho(r,Omega) on all space and angular grid points < this value, the minimizer will consider that a minimum is reach, by default we use 0.001

*`translate_solute_to_center` If true will add Lx/2 (resp Ly/2, Lz/2) to the x(resp y, z) coordinate of the solute read in solute.in WARNING: as described if your solute is not in 0 0 0 this will not translate it at the center of the box.
*`rdfmaxrange`is the value in angstrom until which you want the radial distribution functions of the solute wrt to te solvent to be computed
*`nbinsrdf`  # the radial distribution function is calculated by histograms. nbinsrdf is the number of bins. usually you should have around 4-5 bins per angstrom

*`hydrophobicity` if T you will include long range correction designed to get some hydrophobic behaviour for your solvent, the way it is done is specified in the next keyword
*`treatment_of_hydro` the way you want to compute this hydrophobic correction: VdW use a correction which is based on the Hard Sphere Bridge Correction so if you want to use it you MUST set T to hard_sphere_fluid AND set F to bridge_hard_sphere, the theory is described in jeanmairet et al., THE JOURNAL OF CHEMICAL PHYSICS 139, 154101 (2013). C use a correction that CANNOT BE USE with an HS bridge and it is base on Varrilly et al., THE JOURNAL OF CHEMICAL PHYSICS 134, 074109 (2011) 
*`threebody` if T you will include a correction which will enforse the tetrahedrality of the solvent, this necessary to get good rdf for solvation of ions in water, this is based on the Molinero-Stillinger model described in Molinero et al., J. Phys. Chem. B 2009, 113, 4008–4016. there is currenty two way to add this correction.
*`F3B_old` if T you use a slow way to add the 3 body correction, you should not do that
*`F3B_new`if T you use a faster way to add the 3 body correction, you should always use that if you want to include such a correction.
The parameters of this correction are stored in solute.in, please see the solute.in description.
*`Rc` please ignore this parameter this should be fixed as fast as possible by Guillaume et Maximilien

### solute.in




### solvent.in
