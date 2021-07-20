
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
  - [About Project Source Code](#about-project-source-code)
- [Numerical Method and Algorithms](#numerical-method-and-algorithms)
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
- [Compling](#compling)
    - [Make:](#make)
- [Notes](#notes)

**********

**********

# Introduction
## Origin
&emsp;This project was a coursework for *the Fundamental of Computational Mechanics* of THU in spring 2021. The course contains both CFD and CSD, and lecturers said I could complete a FSI program in place of both final assignments, as a result this program was born. (It was said FSI projects had seldom benn seen in the course apart from one, but I didn't actually get an A for this second record <kbd><del>lmao</del></kbd> )        
## Basic Purpose
&emsp;FSI_T1 is basically a c++ (with some template) library providing multiple classes, which enables you to assemble from them a costom
<kbd style="color:white;background-color:#3F3FBF;font-family:arial;">**[Fluid](#fluid)**</kbd>,&nbsp;
<kbd style="color:white;background-color:#ff9900;font-family:arial;">**[Structural](#structural)**</kbd> 
&nbsp;or&nbsp;
<kbd style="color:white;background-color:#00cc99;font-family:arial;">**[FSI](#fsi)**</kbd>. Fluid part currently supports only Euler equation, which represents adiabatic and invisid compressible ideal gas dynamics. Solid part currently supports linear elastic mechanics, including statics and modal-truncation method dynamics. As the library 



## About Project Source Code
&emsp;All the actual source code files are put in **root** of project directroy, all as .hpp or .h files. The cpp files are meant for case building, containing various case sets.

<br/><br/>

# Numerical Method and Algorithms

## Fluid Methods

## Structural Methods

## FSI Methods

<br/><br/>




# A Brief Guide To Using This Code

## Implementation Framework

## Classes


### Fluid

### Structural

### FSI

## Input File Convention

<br/><br/>

# Case Building

<br/><br/>

# Compling 



### Make:

    make main.exe #for debug
    make mainR.exe #for release
    make mainSG.exe #for mainSG debug
    make mainSGR.exe #for mainSG release

<br/><br/>

# Notes
&emsp;It was only after I STARTED this document when I find my English being rather amateur, even in my own domains... I even searched for a proper antosym for 'proficient' for use in the previous sentence... So please forgive me if any expression in my text seems strange or erroneous.







