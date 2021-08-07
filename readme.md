<!-- markdownlint-disable MD033 MD041-->
<p style="font-weight:bold;font-size:200%;text-align:center;">
FSI_T1
</p>

<p style="font-weight:bold;font-size:150%;text-align:center;">
Fluid-Structure Interaction Practice Homework
</p>

<br><br>

<p style="text-align:center;font-weight:bold;">
by harryzhou2000 @
</p>

<p style="color:#743481;background-color:#FFFFFF;margin-left:20px;padding-left:20px;padding-right:20px;padding-top:20px;padding-bottom:20px;text-align:center;font-size:160%;font-family:Times;">School of Aerospace Engineering, Tsinghua University
</p>

<br><br>

Note:

&emsp;if you speak Chinese, my [report for coursework](docs/三维非结构动网格流固耦合算法设计与实验报告.pdf) is recommended

&emsp;if you desire to go through some analysis of specific cases. If you want to see the latex formulae properly, see [readme.pdf](readme.pdf)

&emsp;this document is under work, the contents are incomplete.

- [Introduction](#introduction)
  - [Origin](#origin)
  - [Basic Purpose](#basic-purpose)
  - [About Project Source Code and Dependencies](#about-project-source-code-and-dependencies)
    - [What did I code?](#what-did-i-code)
    - [Dependencies](#dependencies)
- [Numerical Methods and Algorithms](#numerical-methods-and-algorithms)
  - [Fluid Methods](#fluid-methods)
  - [Structural Methods](#structural-methods)
  - [FSI Methods](#fsi-methods)
- [A Brief Guide To Using This Code](#a-brief-guide-to-using-this-code)
  - [Implementation Framework](#implementation-framework)
    - [Unstructured Meshes](#unstructured-meshes)
    - [FEM Solver](#fem-solver)
    - [FVM Solver](#fvm-solver)
    - [Coupling Solver](#coupling-solver)
  - [Classes](#classes)
    - [Fluid](#fluid)
    - [Structural](#structural)
    - [FSI](#fsi)
  - [I/O File Convention](#io-file-convention)
- [Case Building](#case-building)
- [Compiling](#compiling)
  - [Make](#make)
- [Notes](#notes)
- [References](#references)

**********

**********

# Introduction

## Origin

&emsp;This project was a coursework for *the Fundamental of Computational Mechanics* of THU in spring 2021. The course contains both CFD and CSD, and lecturers said I could complete a FSI program in place of both final assignments, as a result this program was born. (It was said FSI projects had seldom been seen in the course apart from one, but I didn't actually get an A for this second record <kbd><del>lmao</del></kbd> )

## Basic Purpose

&emsp;FSI_T1 is basically a c++ (with some template) library providing multiple classes, which enables you to assemble from them a custom
<kbd style="color:white;background-color:#1F5FCF;font-family:arial;">**[Fluid](#fluid)**</kbd>,&nbsp;
<kbd style="color:white;background-color:#ff9900;font-family:arial;">**[Structural](#structural)**</kbd>
&nbsp;or&nbsp;
<kbd style="color:white;background-color:#00cc99;font-family:arial;">**[FSI](#fsi)**</kbd> computation case. Fluid part currently supports only Euler equation, which represents adiabatic and inviscid compressible ideal gas dynamics. Solid part currently supports linear elastic mechanics, including statics and modal-truncation method dynamics. FSI part, correspondingly, is basically meant for moderate structural displacement, typically aeroelastic cases.

&emsp;The library currently does not have a interface supporting script input or interactive case setting, which means one must construct each single case inside a c++ calling of the relevant classes.

## About Project Source Code and Dependencies

### What did I code?

&emsp;All the actual source code files are put in **root** of project directory, all as .hpp or .h files. The cpp files are meant for case building, containing various case sets, and they can be viewed as some kind of 'static scripts' for the program.

&emsp;The draw-back of not having a script interface and using .cpp as 'script' is that each time you have a new case to compute a re-compiling is needed, which could be rather time consuming.

### Dependencies

&emsp;This project depends on Eigen, which completely is a c++ template library, with no binary file linking. For convenience, I just threw the entire Eigen source code in ./include/, so you won't have to install it. This project also needs c++ standard library and STL to work, which can be handled automatically by compilers mostly.

<br/><br/>

# Numerical Methods and Algorithms

## Fluid Methods

&emsp;To comply with significant movement of boundaries caused by structural displacement, the method of arbitrary Euler-Lagrangian (ALE) description of fluid is adopted. The major differences between ALE and Euler descriptions is that ALE adds some items to Euler relating to the mesh speed[[1]](#ref1).

&emsp;Consider Euler equation for gas dynamics in conservative form:

$$
\frac{\partial{U}}{\partial{t}}+\frac{\partial{F}}{\partial{x}}
+\frac{\partial{G}}{\partial{y}}+\frac{\partial{H}}{\partial{z}}=0\\
U=
\begin{bmatrix}
\rho \\ \rho u \\ \rho v \\ \rho w \\ E  
\end{bmatrix},
F=
\begin{bmatrix}
\rho u\\ \rho u^2 + p \\ \rho uv \\ \rho uw \\ u (E+p)  
\end{bmatrix},
G=
\begin{bmatrix}
\rho v\\ \rho uv \\ \rho v^2 + p \\ \rho vw \\ v (E+p)  
\end{bmatrix},
H=
\begin{bmatrix}
\rho w\\ \rho uw \\ \rho vw \\ \rho w^2 + p \\ w (E+p)
\end{bmatrix}
$$
&emsp;Thus in a ALE control volume:
$$
\frac{\partial}{\partial t}{\int_{\Omega(t)}{UdV}}
+\int_{\partial{\Omega(t)}}{U(u_{\Omega x}n_x+u_{\Omega y}n_y+u_{\Omega z}n_z)d\Gamma}
+\int_{\partial{\Omega(t)}}{(Fn_x+Gn_y+Hn_z)d\Gamma}=0
$$
&emsp;Where $u_{\Omega i}$ are the speed components of the C.V.'s boundaries. The corresponding bold symbols are for vectors in xyz space.

&emsp;Annotating:

$$
\boldsymbol{u^*}=\boldsymbol{u}-\boldsymbol{u_\Omega}
$$

&emsp;Thus:

$$
\frac{\partial}{
  \partial t}{\int_{\Omega(t)}
  \begin{bmatrix}
    \rho \\ \rho u \\ \rho v \\ \rho w \\ E
  \end{bmatrix}
  dV
}
+\int_{\partial{\Omega(t)}}
{
  \left(
    \begin{bmatrix}
    \rho u^*\\ \rho u^*u^* + p \\ \rho u^*v^* \\ \rho u^*w^* \\ u^* (E+p)  
    \end{bmatrix}
    +
    \begin{bmatrix}
    \rho v^*\\ \rho u^*v^* \\ \rho v^*v^* + p \\ \rho v^*w^* \\ v^* (E+p)  
    \end{bmatrix}
    +
    \begin{bmatrix}
    \rho w^*\\ \rho u^*w^* \\ \rho v^*w^* \\ \rho w^*w^* + p \\ w^* (E+p)
    \end{bmatrix}
  \right)d\Gamma
}
\\+\int_{\partial{\Omega(t)}}
{
  \rho \left(
     u_x^* n_x + u_y^* n_y + u_z^* n_z
  \right)
  \begin{bmatrix}
    0 \\ u_{\Omega x} \\ u_{\Omega y} \\ u_{\Omega z} \\ 0
  \end{bmatrix}
}d\Omega
=0
$$

&emsp;Thus if using $\boldsymbol{u^*}$ as speed, the flux could be approximated with a common approximate Riemann solver when the approximate field is not $C^0$. Other items could be handled rather simply.

&emsp;Finite volume spacial discretion is adopted(FVM). This program uses Roe's approximate Riemann solver[[2]](#ref2), which is modified with entropy fix of Harten-Yee[[3]](#ref3) for numerical flux. This program conducts 2nd-order reconstruction with Barth-Jesperson limiter[[4]](#ref4).

&emsp;The boundary conditions are built basically in the form of virtual cells. Slip-wall boundary and non-reflecting boundary are implemented. A pressure inlet/outlet boundary is also implemented , but its robustness is yet to be improved.

&emsp;Time discretion in fluid part is a implicit Euler scheme.

&emsp;The mesh considering FVM is completely composed of tetrahedra.

## Structural Methods

&emsp;Solid part is conducted with classic linear finite element method(FEM). Using symbols form a text book[[5]](#ref5), a linear (small displacement) elastic dynamic problem is discretized as:

$$
M\ddot{a}+C\dot{a}+Ka=F
$$

&emsp;Where a is the vector of nodal DOFs, defined as nodal displacement components. The matrices above are distributed in the finite elements as below:

$$
K^e=\int_{\Omega e}{B^{eT}DB^{e}dV},\ \ \ \;
M^e=\int_{\Omega e}{N^{eT}\rho N^{e}dV},\ \ \ \;\\
F^e=\int_{\Omega e}{fN^edV}+\int_{\partial \Omega e}{pN^ed\Gamma}
$$

&emsp;Where $N^e$ is a row vector of shape functions(interpolating functions)on the nodes in the finite element, and $B^e$ is a matrix representing linear contributions of the nodal DOFs to the strain field(which has 6 rows if using symmetric 3-D strain tensor, and DOF columns corresponding to DOFs). $\rho$ is just density field, and $f$ , $p$ , are external forces of volume and surface.

&emsp;In this project, currently only 4-node tetrahedra elements are implemented, which means both integrals in faces or volumes only need one interpolation point. If higher order elements are to be used, the interpolation method should be accurate enough.

&emsp;A static solver can be found in the program, with LDL decomposition method or PCG method. For FSI purpose, a modal-truncation dynamic solver is implemented. Ignoring matrix $C$ at first, then using a eigen problem:

$$
K\Phi = M\Phi \Lambda
$$

&emsp;While demanding that:

$$
\Phi ^TM\Phi = I
$$

&emsp;Therefore:

$$
\Phi ^TK\Phi = \Lambda
$$

&emsp;Is diagonal.

&emsp;The dynamic system becomes:

$$
\ddot{x}+\Phi ^TC\Phi \dot{x} +\Lambda x=G
$$

&emsp;Where:

$$
G=\Phi ^TF,\ \ \ \ \Phi ^T x = a
$$

&emsp;If the damping part $\Phi ^TC\Phi$ is also diagonal, which is assumed in most CSD computations, the problem is decoupled for each DOF in eigenspace as:

$$
\ddot{x_k}+d_k\dot{x_k}+\lambda _kx_k=g_k
$$

&emsp;Where $d_k$ is a modal damper in the eigenspace, and $\sqrt{\lambda_i}$ are frequencies of each mode.

&emsp;Therefore the discretized problem becomes a set of ODEs that can be separately solved. When the major concern is on vibration rather than delicate propagation of elastic waves, modes with higher frequencies contribute little to major features of the system, so instead of obtaining the whole eigentransformation, modes with n-smallest eigenvalues are solved. As a result, numerical methods solving modes with smallest eigenvalues can be applied here. This project uses the simplest one, the inverse power method.

## FSI Methods

&emsp;The FSI in this project is based on a simple explicit coupling method. The external force for structural model is obtained from the flow in the last coupling time step, and the flow obtains the moving mesh and interface from the structural movement in the last couping time step. If a general implicit method is to be expressed, the condition connecting both problems is the coordinated interface movement and passage of force on the interface, which defines where the overall system should converge. As this project only considers transient problem, only using a explicit method is acceptable.

&emsp;Due to some problems in my meshing technique, the fluid and structural meshes don't adapt perfectly on the interface, but their mesh densities are similar. So when passing force and displacement, this project conducts simple interpolations between the surfaces(with KNN search). It should be noted that this is a common problem in more complex meshing(such as rotating meshes or overlapping grids), so there must be some better and more robust techniques to perform the interpolation than this project's implementation.

&emsp;The last problem in FSI is to move the fluid mesh correctly. There are many interpolating methods that moves the mesh well, but this project uses a much simpler but rather time-consuming way. If you view the fluid part as a solid body, and apply displacement boundary conditions on all its outer faces, when the interface moves, the fluid mesh moves correspondingly. When the interface is attached, the movements inside the mesh is also almost continuous. The only problem that arises is that 'stress' distribution in the fluid mesh could be singular around some points of the interface, and is mostly concentrated near the interface. Under the assumption of linear structural response, the elements near the interface are likely to be overlapped or distorted. A simple way of solving this is to set higher value of modulus near the interface, forcing the strain to distribute farther from the interface. Also, setting higher shearing modulus could help abating the distortion.

&emsp;To set the modulus distribution automatically, a distance field could be applied. However, in this project, as the elements are actually smaller near the interface, the volumes of interface are used to set the modulus.

<br/><br/>

# A Brief Guide To Using This Code

## Implementation Framework

### Unstructured Meshes

&emsp;A good representation of the topology of an unstructured mesh, would be a diagram representing hierarchical connectivity between layers of mesh elements, which is frequently referred to as a *Hasse diagram*. For example, layers could be **vertices** $\rArr$ **edges** $\rArr$ **faces** $\rArr$ **volumes** for a 3-D mesh, while in some certain problems, some intermediary layers could be omitted. The diagram is essentially a DAG whose nodes represent mesh elements. Therefore, the corresponding data structure would be basically a sparse adjacency matrix, or separate matrices between levels of layers. With this kind of information, one can easily find out, for example, the set of vertices sharing edge, face or element with a known vertex, or the set of volumes sharing common faces or nodes with a known volume. These queries naturally take $O(1)$ time with the help of the Hasse diagram.

&emsp;Normally, the geometry of a mesh is represented as known coordinates of vertices.

&emsp;When I started this project, it was actually based on some previous work with some rather dumb coding. From the view of a Hasse diagram, my unstructured mesh is only representing relations between volumes and vertices, which is alright with FEM coding, for neither volume and surface integrals nor nodal average in FEM need any more information. In 2nd order FEM, neighboring cells are need, for DOFs are stored with volumes and their adjacency is the same as cellular adjacency. I did not upgrade my data structure to represent one more layer (the face layer), but only added a neighbor searching procedure to the mesh loading method. Apparently, this kind of technique prevents direct extension to higher-order FVM with large stencils. Also, the adjacency representation is not generic, which only supports tetrahedral volumes(cells) with only 4 vertices, which hinders further implementation of higher-order nodal FEM.

&emsp;In namespace **MeshBasic**, class **gridNeighbour** is defined, which records information about neighboring cells for a single cell, including some precalculated geometric values. Mean, class **TetraNodes** is also defined here to store actual volume to vertices info for a single cell, along with some precalculated geometric information. In the same namespace, I defined class **TetraMesh**, which is a general mesh holder. TetraMesh is able to load meshes and precalculate necessary topological and geometric information (which is primarily for FVM).

&emsp;I've established some understandings of unstructured mesh, those are recorded below and irrelevant with the current project.

&emsp;As I truly know little about unstructured mesh management and related data structures, I can only make some assumptions on a nearly 'perfect' implementation of an unstructured mesh framework. This framework should be able to represent a general mesh topology, irrelevant with mesh type, node distribution and dimension. More importantly, the framework should provide methods of obtaining adjacent mesh elements with constant time complexity (with degrees of the graph limited). Also, the frame work should allow users to attach arbitrary form of data onto mesh elements, to represent scalar, vector or tensor field within the mesh geometry, on vertices, faces or volumes. Meanwhile, the framework should enable users to abstract these fields into global vectors and tensors, in cooperation with global discrete operators like sparse matrices. Of course, process-scale parallel is a must, most probably implemented with MPI framework. Thread parallelism and CUDA support are also needed, but message distributing and communication are still the most complicated part of the program.

### FEM Solver

&emsp;In conventional finite-element method, or a common Galerkin method, the trial functions, being identical with the bases of the discrete solution, are considered to be $C^0$. Therefore the piecewise defined polynomial bases should maintain $C^0$ continuity on the interfaces of volumes (take 3-D for example). As a result, it is better to consider the discrete DOFs to be set on the vertices, and the bases are interpolation functions that satisfy a kronecker-delta property over the vertices. For this certain project, each volume is a tetrahedron with 4 vertices, and you can easily transform it linearly into a corner of a cartesian box, and using 3 cartesian axes you can easily define first-order bases functions for each vertex. Generally, the basis functions are defined in a normalized coordinate system as polynomials. As the transformation between the normalized space and the geometric space is generally not a linear mapping (basically as curved elements), generally the same set of basis functions are used to interpolate the mapping. Fundamentally, the actual bases are fractions rather than polynomials, but as flat-faced and near flat-faced elements are the majority, the orders of polynomial-based numerical integral are mostly decided with the situation of a linear spacial mapping.

&emsp;Sadly, in this program you may find no actual numerical integral process, as the variable to integrate are derived form first derivatives of the bases, and they are actually constant in the volumes for linear bases. The only extra value to calculate for each volume is the determinant of the first partial derivatives (or Jacobian matrix) of the spacial mapping (which is also constant in the volume). The determinant is precisely the volume of the tetrahedron.

&emsp;When discussing a linear elastic static problem, all should be calculated are the discrete linear operator (or stiffness matrix) and the load vector. Using some linear algebra techniques, the program first integrates a local matrix and a local load vector for each volume, and adds them to the global matrix and global load vector. This program applies triplet structure to the sparse matrix, by adding random entries first and doing a sort on the indices in the end, it needs $O(N logN)$ time for assembly (mainly for sorting). Certainly, taking advantage of the (at least temporarily) static feature of mesh, a preallocation procedure could be added and one can easily build the matrix in CSR form. The CSR preallocation-and-fill paradigm takes $O(N log(Degree))$ time, where Degree denotes an average number of non-zeros in a row. This latter kind of implementation cannot be found in the project for it's a few times more complicated.

&emsp;There are some other problems concerning boundary conditions, which are small modifications in the procedures above (although actually the load vector's non-zero entries are mainly caused by the boundary conditions) concerning reducing DOFs or boundary integration.

&emsp;When the matrix and vector are produced, the only matters are to solve them. For a time-dependent problem, you need to derive a proper time discretion scheme; for a non-linear problem, you need to update the stiffness matrix at certain times. Nonetheless, the core problem here is to solve the linear system. I just threw them to Eigen's internal template implementation (While being a C++ template library, Eigen is also a great C++ interface and it supports wrapping of various external solvers). For a elliptic problem here with Galerkin discretion, matrices are mostly positive-definite and symmetric, so PCG would be rather favorable. When the problem is not too large, direct methods like LDU decomposition could also be viable (which is extremely helpful in eigenvalue problems for they require a lot of re-solving the matrix). Concerning the sparse eigenvalue problem (for modal analysis in FSI), I simply applied inverse power method, which seems to have some problems in high-rank mode convergence. Should consider switching to some more advanced techniques or just use a well-proven library.

&emsp;The global vectors are represented as std::vector\<T\>, for I only apply shared-memory parallelism. All the matrix-related data structures and algorithms are in namespace **SparseMat**, and the **SparseMatGeR** class with relevant solving functions inside. **SolidMaterial** namespace and **ElasSet** class defines some general constitutional properties for the elastic problem. In **FEM** namespace class **FemT4Solver** is defined, which is derived from the **TetraMesh** class (which I now think should become a member rather than father... but inheritance means all the complex data in TetraMesh could be written the same way in TetraMesh member functions...). The **FemT4Solver** class, with **TetraMesh**, imports and manages mesh and problem definition, and assembles the stiffness matrix along with load vector, and provide interface to solve and output.

<!--TODO-->

### FVM Solver

### Coupling Solver

## Classes

### Fluid

### Structural

### FSI

## I/O File Convention

<br/><br/>

# Case Building<!--TODO-->

<br/><br/>

# Compiling

&emsp;In theory, this project should be included as a header-only c++ library in your project, which means no actual building is needed. However, as the input/output and intermediate files are pretty tricky to produce, so main.cpp is provided for using. mainbkp.cpp provides other useful cases than in main.cpp. To compile, you simply use **g++** to compile the single main.cpp file, or use [**make**](#make) like instructed below.

&emsp;If you put all the functions in mainbkp.cpp into main.cpp to compile, when using mingw64 g++ in Windows, the compiler could exit abnormally declaring an error on insufficient resources concerning templates. I can't repeat this problem, please tell me in the issues why this could happen if you know.

&emsp;The mainSG.cpp is for a 2D structural gird Euler equation gas dynamics solver, sharing some headers. Its case construction appears much more complex.

<!--TODO-->

## Make

    make main.exe #for debug
    make mainR.exe #for release
    make mainSG.exe #for mainSG debug
    make mainSGR.exe #for mainSG release

<br/><br/>

# Notes

&emsp;It was only after I STARTED this document when I find my English being rather amateur, even in my own domains... I even searched for a proper antonym for 'proficient' for use in the previous sentence... So please forgive me if any expression in my text seems strange or erroneous.

<br/><br/>

# References

<span id="ref1"></span>
[1] P, Le, Tallec, et al. Fluid structure interaction with large structural displacements[J]. *Computer Methods in Applied Mechanics and Engineering*, 2001, 190(24-25):3039-3067.

<span id="ref2"></span>
[2] Roe P L . Approximate Riemann solvers, parameter vectors, and difference schemes[J]. *Journal of Computational Physics*, 1981, 43(2):357-372.

<span id="ref3"></span>
[3] Yee H C . Upwind and symmetric shock-capturing schemes. 1987.

<span id="ref4"></span>
[4] Barth T J ,  Jespersen D C . The design and application of upwind schemes on unstructured meshes[J]. AIAA Aerospace Sciences Meeting, 1989, 0366(13).

<span id="ref5"></span> [5] 王勖成. 有限单元法[M]. 清华大学出版社, 2003.
