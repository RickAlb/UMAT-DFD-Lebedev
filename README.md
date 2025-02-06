# UMAT-DFD-Lebedev

## Purpose
Abaqus-Standard subroutine based on the Discrete Fibers Dispersion (DFD) model in finite strains [(Ehret et al., 2010; Li et al., 2018)](#references) for fibre-reinforced soft biological materials.
Features of the model:

- isotropic ground matrix with neo-Hookean model;
- anisotropic fibre contribution with multiple fiber families and independent in-plane and out-of-plane fiber dispersion;
- exponential fibre strain energy;
- fiber families principal basis {$`\mathrm{M}_i`$, $`\mathrm{M}_{\mathrm{ip}i}`$, $`\mathrm{M}_{\mathrm{op}i}`$} definable relative to:
  * global basis {$`\mathrm{G}_1`$, $`\mathrm{G}_2`$, $`\mathrm{G}_3`$};
  * local basis {$`\mathrm{E}_1`$, $`\mathrm{E}_2`$, $`\mathrm{E}_3`$};
- discrete integration:
  * Lebedev quadrature (built-in);
  * integration points from external file;
- volumetric behavior formulations:
  * incompressible (2D plane stress, 3D hybrid formulation);
  * nearly incompressible (3D hybrid total Lagrange multiplier formulation);
  * nearly incompressible (3D decoupled deviatoric and volumetric strain-energies);
- uniform prestress field $`\sigma_0`$ applicable to all elements in the model (the components of the relative deformation gradient $`\mathrm{F}_0`$ must be provided through a dedicated prestretch text file).

The model is implemented applying the rotation correction to the deformation gradient $`\mathrm{F}`$ when a local orientation {$`\mathrm{E}_1`$, $`\mathrm{E}_2`$, $`\mathrm{E}_3`$} is specified in conjunction with Solid (Continuum) elements. See [Nolan et al. (2022)](#references) for further details.

Includes the subroutine for the Lebedev integration points computation [(Lebedev & Laikov, 1999)](#references).

Special functions subroutines are also included under Shanjie Zhang and Jianming Jin permission [(Zhang & Jin, 1996)](#references).

## Pointers setup
The subroutine uses global arrays to store the integration points of each defined DFD material. The arrays can be read from any thread at any time during the analysis by means of pointers.
To ensure the subroutine compiles correctly, open the environment file `abaqus_v6.env` in the Abaqus installation folder and add the following flags to the Fortran compiler settings:

- Linux platform (GCC compiler): ` -fcray-pointer`;
  * For example `compile_fortran = (fortCmd + " -free" + " -fcray-pointer" + " -c -fPIC -I%I")`
- Windows platform (Intel Fortran Compiler): no action needed.

The subroutine includes the Fortran interface files `SMAAspUserArrays.hdr` and `SMAAspUserUtilities.hdr` required to create and access the pointers (see [2.1.23 Allocatable arrays, Abaqus User Subroutines Reference Guide]()).
For the older versions of Abaqus, e.g. Abaqus 6.14, these files are located in `%AbqInstall%\Abaqus\6.14-1\code\include`.
In Abaqus 2019, these files are located in `%AbqInstall%\SimulationServices\V6R2019x\win_b64\code\include`, but contain an incorrect `#` character. To make them work properly, remove the extra `#` and save.

## Abaqus input material definition
The subroutine `UMAT_DFD_LEB.for` supports up to 100 different DFD materials defined at the same time during the analysis.
To define a DFD material add the following lines to the Abaqus input file:

```
*MATERIAL, NAME=MATNAME
*USER MATERIAL, TYPE=MECHANICAL, CONSTANTS=N
  bulk,    νg,     μ,    c1,    c2,   νf1,    α1,    β1,
    γ1,    a1,    b1,   ...,   ...,   ...,   ...,   ...,
   νfm,    αm,    βm,    γm,    am,    bm,  rint,  nint,
  thrp,   csw,  elem,   ori,    θ1,    θ2,    θ3
```

where `MATNAME` is the material name to be defined according to [Material name and volumetric behavior](#material-name-and-volumetric-behavior) section, and `N` is the total number of parameters defined.
Parameters must be provided 8 per line.
See [*USER MATERIAL, Abaqus Keywords Reference Guide]() for further information.

### Material name and volumetric behavior
The material name `MATNAME` must be max 80 characters long, and must be defined as `DFD_XX-custom_material_name` where `XX` defines the volumetric behavior.
Available options are:

- `IN` for incompressible formulaton;
- `NI` for nearly incompressible hybrid total Lagrange multiplier formulation.
- `HB` nearly incompressible formulation with decoupled deviatoric and volumetric strain-energies. The deviatoric fiber strain-energy (DFD model) is implemented using the full deformation gradient according to [Helfenstein et al. (2010)](#references) and [Nolan et al. (2014)](#references).

### Elements
The DFD model is implemented using the objective fourth-order tangent stiffness tensor relative to the Jaumann rate of the Kirchhoff stress.
This stiffness tensor is intended to work in conjunction with Solid (Continuum) elements (see Table 1.5.3–1, [1.5.3 Stress rates, Abaqus Theory Guide]()), but Structural elements can also be used with some limitations (Membranes only).
Elements are supported in full integration only, and can be selected depending on the volumetric behavior as follows:

- `IN`:
  * 2D plane stress (Continuum) elements `CPS3`, `CPS4`, `CPS6`, `CPS8`;
  * 2D membrane (Structural) elements `M3D3`, `M3D4`, `M3D6`, `M3D8`, `M3D9`;
  * 3D hybrid (Continuum) elements `C3D4H`, `C3D6H`, `C3D8H`, `3D10H`, `C3D15H`, `C3D20H`. Add the option `HYBRID FORMULATION=INCOMPRESSIBLE` to `*USER MATERIAL` to activate the hybrid formulation;
- `NI`:
  * 3D hybrid (Continuum) elements `C3D4H`, `C3D6H`, `C3D8H`, `C3D10H`, `C3D15H`, `C3D20H`. Add the option `HYBRID FORMULATION=TOTAL` to `*USER MATERIAL` to activate the hybrid formulation;
- `HB`:
  * 3D (Continuum) elements `C3D4`, `C3D6`, `C3D8`, `C3D10`, `C3D15`, `C3D20`.

### Material and computation parameters
Parameters must be given in the following order:

| no.  | Parameter | Unit | Description |
| ------------- | ------------- | ------------- | ------------- |
| $`1`$ | `bulk` | Pa, kPa, ... | Bulk modulus (ignored in incompressible formulation, `IN`). |
| $`2`$ | `νg` | - | Ground matrix volume fraction. |
| $`3`$ | `μ` | Pa, kPa, ... | Ground matrix shear modulus. |
| $`4`$ | `c1` | Pa, kPa, ... | Fibers stiffness-like parameter. |
| $`5`$ | `c2` | - | Fibers stiffeneing parameter. |
| $`6`$ | `νf1` | - | Fiber family's volume fraction. |
| $`7`$ | `α1` | deg | Azimuth angle of the fiber family's principal basis {$`\mathrm{M}_1`$, $`\mathrm{M}_{\mathrm{ip}1}`$, $`\mathrm{M}_{\mathrm{op}1}`$} from the global basis vector $`\mathrm{G}_1`$ or the local basis vector $`\mathrm{E}_1`$ if active. |
| $`8`$ | `β1` | deg | Elevation angle of the fiber family's principal basis {$`\mathrm{M}_1`$, $`\mathrm{M}_{\mathrm{ip}1}`$, $`\mathrm{M}_{\mathrm{op}1}`$} from the global basis plane {$`\mathrm{G}_1`$, $`\mathrm{G}_2`$} or the local basis plane {$`\mathrm{E}_1`$, $`\mathrm{E}_2`$} if active. |
| $`9`$ | `γ1` | deg | Rotation of the fiber family's principal basis {$`\mathrm{M}_1`$, $`\mathrm{M}_{\mathrm{ip}1}`$, $`\mathrm{M}_{\mathrm{op}1}`$} about $`\mathrm{M}_1`$. |
| $`10`$ | `a1` | - | In-plane von Mises distribution concentration parameter, $`-\infty < a < \infty`$ (ignored if `rint` is `999`). |
| $`11`$ | `b1` | - | Out-of-plane von Mises distribution concentration parameter, $`-\infty < b < \infty`$ (ignored if `rint` is `999`). |
| ... |  |  | Add $`n`$ fiber families by repeating parameters from $`6`$ to $`11`$. |
| $`6n+6`$ | `rint` |  | Integration scheme flag. Set `0` to use implemented Lebedev quadrature, or `999` to read points from external file(s). |
| $`6n+7`$ | `nint` |  | Target number of Lebedev integration points on the unit sphere **per fiber family** (ignored if `rint` is `999`). The analysis will be performed using the set containing the number of points closer to `nint`. Then, only the points on the top hemisphere will be used. |
| $`6n+8`$ | `thr` |  | Threshold value for integration points retention. Points with integrand value lower than `thr` will be discarded. Set `1e-6` for best performance/integration accuracy, `0` to keep all positive integration points, or else for custom computation. **Warning:** discarding too much points with excessively high threshold produces inaccurate results.|
| $`6n+9`$ | `csw` |  | Fibers compression switch flag. Set `1` to exclude compressed fibers, or `0` to include all fibers. |
| $`6n+10`$ | `elem` |  | Element type. Set `1` for Continuum elements, or `2` for Structural elements. |
| $`6n+11`$ | `ori` |  | Type of active orientation. Set `0` if global orientation is used, {$`\mathrm{G}_1`$, $`\mathrm{G}_2`$, $`\mathrm{G}_3`$} basis, or `1` if local orientation is used, {$`\mathrm{E}_1`$, $`\mathrm{E}_2`$, $`\mathrm{E}_3`$} basis. |
| $`6n+12`$ | `θ1` | deg | Rotation of {$`\mathrm{E}_1`$, $`\mathrm{E}_2`$, $`\mathrm{E}_3`$} about $`\mathrm{G}_1`$ (ignored if `ori` is `0`). |
| $`6n+13`$ | `θ2` | deg | Rotation of {$`\mathrm{E}_1`$, $`\mathrm{E}_2`$, $`\mathrm{E}_3`$} about $`\mathrm{G}_2`$ (ignored if `ori` is `0`). |
| $`6n+14`$ | `θ3` | deg | Rotation of {$`\mathrm{E}_1`$, $`\mathrm{E}_2`$, $`\mathrm{E}_3`$} about $`\mathrm{G}_3`$ (ignored if `ori` is `0`). |

Notes:
- Stiffness parameters must be expressed in the unit consistent with the force-area units used in the model.
- The material must include at least one fiber family ($`n=1`$);
- The sum of all the fiber families volume fractions `νfi` should be equal to 1;
- Parameters `elem` and `ori` are essential for addressing the rotation correction of the deformation gradient $`\mathrm{F}`$ [(Nolan et al., 2022)](#references);
- The rotation angles `θ1`, `θ2`, `θ3` of the local basis {$`\mathrm{E}_1`$, $`\mathrm{E}_2`$, $`\mathrm{E}_3`$} are used to print the  components of the deformation gradient $`\mathrm{F}`$ and the Cauchy stress tensor $`\sigma`$ relative to the global basis {$`\mathrm{G}_1`$, $`\mathrm{G}_2`$, $`\mathrm{G}_3`$} on a `.chy` file in the analysis' folder. The target element and the relative Gauss point number must be specified in the `kpars_umat.for` source file. 

### Reading integration points from file
When the integration scheme flag is `rint` is `999`, the integration points must be given on differet files for each fiber family using the material name `MATNAME` as file name, and the progrssive file extension `.intp[i]` to identify the ith family.
For example, a material named `DFD_IN-example_material` with `N` fiber families must be provided with the following files:
- `%AnalysisFolderPath%\DFD_IN-test_material.intp1` for the 1st fiber family;
- `%AnalysisFolderPath%\DFD_IN-test_material.intp2` for the 2nd fiber family;
- ...
- `%AnalysisFolderPath%\DFD_IN-test_material.intpN` for the Nth fiber family;

A reference `.intp[i]` file is given along with the input files in SET LINK TO DIRECTORY 

## Prestress definition (optional)
Prestress is supported **ONLY** when the following two conditions are met:
1. prestress is uniformly applied to all the elements in the model;
2. only one DFD material is defined.

Prestress applied to a subset of elements or applied in presence of more than one DFD material is **NOT** supported.

To define a uniform prestress field add the following lines to the Abaqus input file:

```
*INITIAL CONDITIONS, TYPE=STRESS
ALL_ELEMENTS_SET, σ11, σ22, σ12,
```

where `σ11`, `σ22`, `σ12` are the components of the initial stress field $`\sigma_0`$ for a plane stress problem.
Refer to [*INITIAL CONDITIONS, Abaqus Keywords Reference Guide]() for further information about prestress in general 3D problems.

To allow the correct computation of the tangent stiffness matrix, the initial deformation gradient $`\mathrm{F}_0`$ associated to the prestress field $`\sigma_0(\mathrm{F}_0)`$ must be provided through a text file `INPUTFILENAME.F0` included in the analysis' folder, where `INPUTFILENAME` must be the same as the analysis input file.
The components must be formatted as exponential with 16 positions and 6 decimal digits (fortran format `E16.6`), and organized as a matrix on three rows and columns. For example, for a uniaxial prestress field along $`\mathrm{G}_2`$ (or $`\mathrm{E}_2`$), $`\sigma_{11}=\sigma_{12}=0`$, $`\sigma_{22}>=0`$, $`\mathrm{F}_0`$ must be written as follows:

```
    8.560503E-01    0.000000E+00    0.000000E+00
    0.000000E+00    1.169145E+00    0.000000E+00
    0.000000E+00    0.000000E+00    9.991537E-01
```

## Launch analysis
Use the following command line to launch the job from terminal:

```
abaqus job=input_file_name user=main_sub cpus=j
```

where `main_sub` is the main subroutine `main_sub.for` including `UMAT_DFD_LEB.for` and all the subroutines required for the analysis. The subroutine supports parallel computing. Replace `j` with the desired number of cores (default is `1` by omitting `cpus` option).
See [3.2.2 Abaqus/Standard execution, Abaqus Analysis User's Guide]() for further details.

## References

Ehret, A. E., M. Itskov, and H. Schmid. “Numerical Integration on the Sphere and Its Effect on the Material Symmetry of Constitutive Equations-A Comparative Study.” International Journal for Numerical Methods in Engineering 81, no. 2 (January 8, 2010): 189–206. https://doi.org/10.1002/nme.2688.

Helfenstein, J., M. Jabareen, E. Mazza, and S. Govindjee. “On Non-Physical Response in Models for Fiber-Reinforced Hyperelastic Materials.” International Journal of Solids and Structures 47, no. 16 (August 1, 2010): 2056–61. https://doi.org/10.1016/j.ijsolstr.2010.04.005.

Lebedev, Vyacheslav Ivanovich, and Dmitri Laikov. “A Quadrature Formula for the Sphere of the 131st Algebraic Order of Accuracy.” In Doklady Mathematics, 59:477–81. Pleiades Publishing, Ltd., 1999.

Li, Kewei, Ray W. Ogden, and Gerhard A. Holzapfel. “A Discrete Fibre Dispersion Method for Excluding Fibres under Compression in the Modelling of Fibrous Tissues.” Journal of The Royal Society Interface 15, no. 138 (January 31, 2018): 20170766. https://doi.org/10.1098/rsif.2017.0766.

Nolan, D.R., A.L. Gower, M. Destrade, R.W. Ogden, and J.P. McGarry. “A Robust Anisotropic Hyperelastic Formulation for the Modelling of Soft Tissue.” Journal of the Mechanical Behavior of Biomedical Materials 39 (November 2014): 48–60. https://doi.org/10.1016/j.jmbbm.2014.06.016.

Nolan, D.R., C. Lally, and J.P. McGarry. “Understanding the Deformation Gradient in Abaqus and Key Guidelines for Anisotropic Hyperelastic User Material Subroutines (UMATs).” Journal of the Mechanical Behavior of Biomedical Materials 126 (February 1, 2022): 104940. https://doi.org/10.1016/j.jmbbm.2021.104940.

Zhang, Shanjie, and Jianming Jin. Computation of Special Functions. Wiley, 1996., ISBN 0-471-11963-6, LC QA351.C45.


