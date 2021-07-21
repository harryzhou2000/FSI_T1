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
  - [Classes](#classes)
    - [Fluid](#fluid)
    - [Structural](#structural)
    - [FSI](#fsi)
  - [Input File Convention](#input-file-convention)
- [Case Building](#case-building)
- [Compiling](#compiling)
  - [Make](#make)
- [Notes](#notes)
- [References](#references)

**********

**********

# Introduction

&emsp;Note: if you speak Chinese, my [report for coursework](docs/三维非结构动网格流固耦合算法设计与实验报告.pdf) is recommended if you desire to go through some analysis of specific cases.

## Origin

&emsp;This project was a coursework for *the Fundamental of Computational Mechanics* of THU in spring 2021. The course contains both CFD and CSD, and lecturers said I could complete a FSI program in place of both final assignments, as a result this program was born. (It was said FSI projects had seldom been seen in the course apart from one, but I didn't actually get an A for this second record <kbd><del>lmao</del></kbd> )

## Basic Purpose

&emsp;FSI_T1 is basically a c++ (with some template) library providing multiple classes, which enables you to assemble from them a custom
<kbd style="color:white;background-color:#1F2F6F;font-family:arial;">**[Fluid](#fluid)**</kbd>,&nbsp;
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
+\int_{\partial{\Omega(t)}}
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

&emsp;Finite volume spacial discretion is adopted(FVM). This program uses Roe's approximate Riemann solver[[2]](#ref2), which is modified with entropy fix of Harten-Yee[[3]](#ref3) for numerical flux. This program conducts 2nd-order
reconstruction with Barth-Jesperson limiter[[4]](#ref4).

&emsp;The boundary conditions are built basically in the form of virtual cells. Slip-wall boundary and non-reflecting boundary are implemented. A pressure inlet/outlet boundary is also implemented , but its robustness is yet to be improved.

&emsp;Time discretion in fluid part is a implicit Euler scheme.

&emsp;The mesh considering FVM is completely composed of tetrahedra.

## Structural Methods<!--TODO-->

## FSI Methods<!--TODO-->

<br/><br/>

# A Brief Guide To Using This Code<!--TODO-->

## Implementation Framework

## Classes

### Fluid

### Structural

### FSI

## Input File Convention

<br/><br/>

# Case Building<!--TODO-->

<br/><br/>

# Compiling<!--TODO-->

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
