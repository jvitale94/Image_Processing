# Image_Processing
An image processing library with custom methods and filters. I took a computer vision course at the Aquincum Institue of Technology in the Fall of 2015, and this is the code I wrote for that course, with some extra messing around thrown in there. I only implemented the methods in the R2Image.cpp file in the src directory, from line 250 to line 1572. The most important method in this file to the blendOtherImageHomography method beginning on line 1243, as it uses most pervious methods to blen two offset images together. I have included a sample output of this method in the output folder.

Everything else in the skeleton directory was given to me by my professor at the Aquincum Institue of Technology, Gergely Vass. All credit for code other than what I have already specified goes to him. I have also included the README that he gave my class below, which gives instructions of how to compile the code. To run the code, go to the skeleton directory and paste the following command: 
`src/imgpro input/bridge01.jpg output/test.jpg -matchHomography input/bridge02.jpg ` 

This is a sample, as you can call many other methods. Use the -help commands to see the other methods you can call.

##Gergely Vass's README

This directory contains skeleton code as a starting point for Computer Vision. 


FILE STRUCTURE
==============

There are several files, but you should mainly change src/R2Image.cpp.

src/ - Directory with source code
Makefile - Unix/Mac makefile for building the project with "make". 
imagepro.[vcproj/sln/suo] - Project file for Visual Studio 2005 on Windows
imgpro.cpp - Main program, parses the command line arguments, and calls the appropriate image functions
R2Image.[cpp/h] - Image class with processing functions (this is the only file that you need to edit)
R2Pixel.[cpp/h] - Pixel class 
R2/ - A library of useful 2D geometric primitives
jpeg/ - A library for reading/writing JPEG files
input/ - Contains example input images. 
output/ - Es empty to start -- it will contain the images produced by your program (see below)
runme.bat - a script (for Windows) that you will fill in to demonstrate execution of your program
runme.sh - same as <code>runme.bat, but for Mac OS X

COMPILATION
===========

If you are developing on a Windows machine and have Visual Studio
installed, use the provided project solution file (assn1.sln) in the
src/ directory to build the program. If you are developing on a Mac or
Linux machine, cd into the src/ directory and type "make". In either
case, an executable called imgpro (or imgpro.exe) will be created in
the src/ directory.
