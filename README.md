**FEMFFUSION**: a finite element method code for nuclear reactor modelling  

```
      ___  ___  _   _  ___  ___  _ _  __  _  _  _  _   
     | __|| __|| \_/ || __|| __|| | |/ _|| |/ \| \| |  
     | _| | _| | \_/ || _| | _| | U |\_ \| ( o ) \\ |  
     |_|  |___||_| |_||_|  |_|  |___||__/|_|\_/|_|\_|  
     
                    BASIC VERSION   1.0  
```
[**FEMFFUSION**](https://bitbucket.org/Zonni/femffusion/) is an open source C++ neutronic code that solves the multigroup neutron transport equation using the diffusion approximation or the SPN approximation. The code uses the continuous Galerkin finite element method to be able to deal with any type of geometry and problem dimension (1D, 2D and 3D problems). It works on top of [deal.II](https://www.dealii.org/) library, which provides supporting and advances in the finite element method. 

**FEMFFUSION** main features are:

 * Open source software (released under the terms of the [GNU GPL version 3](http://www.gnu.org/copyleft/gpl.html))
 * Use of the FEM to solve the multigroup neutron diffusion (or SPN) equations.
 * Use of matrix-free technique to maintain reasonable memory demands.
 * Possibility to use and implement a variety of eigenvalue solvers and preconditioners. 
 * Valid for all types of geometries, rectangular, hexagonal, pin-level and unstructured.
 * Possibility to import unstructured grids from [GMsh](http://gmsh.info/).
 * Capacity to solve problems in 1D, 2D and 3D.
 * Solve the direct or the adjoint flux and several eigenpairs, if requested.
 * Time-dependent problems as: control rod movements, noise problems, custom time-dependent cross sections…
 * Frequency domain analysis of noise problems.
 * Output provides the keff, the map of the averaged neutron power per assembly and each of the fluxes. Also standard .vtk files are provided among other output formats. For time dependent problems the  neutron power and the neutron fluxes are provided each time step. 
 * Interact with high-quality open-source libraries: deal.II, PETSc, SLEPc, Sundials...
 * Easy interface with plotting and post-processing tools (Matlab, ParaView, Matplotlib...).
 * Well documented and easy to extend to related problems.


Actually, *FEMFFUSION* can be seen as deal.II specialization for the neutron transport equation that uses large sparse matrix solvers ([SLEPc](http://www.grycap.upv.es/slepc/) among others). That is to say, *FEMFFUSION* builds the matrices ![](https://latex.codecogs.com/gif.download?L) and ![](https://latex.codecogs.com/gif.download?M) that cast the multigroup neutron transport/diffusion equation as a matrix-based eigenvalue problem:

![my equation](https://latex.codecogs.com/gif.download?L%20%5Cphi%20%3D%20%5Cfrac%7B1%7D%7Bk_%5Ctext%7Beff%7D%7D%20M%20%5Cphi)

These matrices are expected to be sparse, as they are the result of the discretization of the differential transport operator using the finite element method, over a certain spatial grid either. Said matrices are thus built in [PETSc](http://www.mcs.anl.gov/petsc/) format, so they can either be passed to a solver (default is [SLEPc](http://www.grycap.upv.es/slepc/), whose algorithms and parameters may be chosen at run-time). 

*FEMFFUSION* also provides a second glue layer that links the output of the linear/eigen-solver to the input of a post-processing tool ([ParaView](http://www.paraview.org/) ). The effective multiplication factor, keff, is shown in the output file along with the fluxes and power distribution.
 


# Install

*FEMFFUSION* works under GNU/Linux. For Debian based distributions an installing script is provided.
First, install some important tools and libraries by typing in the terminal:

```
sudo apt-get install make autoconf automake gcc g++ git findutils cmake
sudo apt-get install valgrind libgsl-dev petsc-dev slepc-dev
sudo apt-get install gmsh gnuplot paraview 
```

Then, clone the *FEMFFUSION* repository where you want it, dowload other libraries, compile and check. This will take some time due to compilation of pestc, slecp and deal.II, around 1 hour. It requires Internet connection.

```
git clone https://Zonni@bitbucket.org/Zonni/femffusion.git
cd femffusion
./install.sh
```

If you have any problem, you can contact <anvifer2@upv.es>.

# Examples 

After the successful compilation of the code one recommended step is to run all the examples by:

```
./run_examples.sh
```

It consists of a set of nuclear reactors of 

 1. **1D**; Example of a homogeneous slab of 2cm solved with vaccum, zero current and zero flux boundary conditions. A pedagogical example to see the improvement between diffusion, SP3 and SP5 and transport methods.

| ![1D](doc/figures/1.png){width=50}   |
|:----------:|
| 1st Mode   |  


 2. **2D_BIBLIS**; Classic neutron difusion benchmark with a characteristic chess board pattern. Solved with SP1, SP3 and SP5 aproximations.

  
| ![FirstMode](doc/figures/BIBlISMPc1.png)  |    ![SecondMode](doc/figures/BIBlISMPc2.png)|
|:----------:|:-------------:|
| 1st Mode |  2nd Mode  |
| ![ThirdMode](doc/figures/BIBlISMPc3.png)  |    ![FourthMode](doc/figures/BIBlISMPc4.png)|
| 3rd Mode |  4th Mode  |
 
 3. **2D_C5G7**: The NEA 2D benchmark on deterministic transport
calculations without spatial homogenisation. The description can be found [here](https://www.oecd-nea.org/science/docs/2003/nsc-doc2003-16.pdf).

| ![FirstMode](doc/figures/c5g5_power_1-1.png)  |    ![SecondMode](doc/figures/c5g5_power_2-1.png)|
|:----------:|:-------------:|
| 1st Mode |  2nd Mode  |
| ![ThirdMode](doc/figures/c5g5_power_3-1.png)  |    ![FourthMode](doc/figures/c5g5_power_3-1.png)|
| 3rd Mode |  4th Mode  |

 4. **2D_VVER440**: An hexagonal two dimensional benchmark.

 5. **3D_IAEA**: A classic 3D benchamark with rectangular geometry. The description can be found [here](https://engineering.purdue.edu/PARCS/Code/TestSuite/CalculationMode/StandAloneMode/Eigenvalue/IAEA3DPWR).
 
 6. **3D_Langenbuch**: Other classic 3D benchamark.

 7. **3D_VVER440**: An hexagonal 3D benchmark.



# Further information


Repository: <https://bitbucket.org/Zonni/femffusion> 
Authors: Antoni Vidal-Ferràndiz, Amanda Carreño, Damian Ginestar and Gumersindo Verdú.
Universitat Politècnica de València
Contact: <anvifer2@upv.es>


*FEMFFUSION* is licensed under [GNU GPL version 3](http://www.gnu.org/copyleft/gpl.html).
*FEMFFUSION* is  [free software](https://www.gnu.org/philosophy/free-sw.html): you are free to change and redistribute it.
There is NO WARRANTY.

