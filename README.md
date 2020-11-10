# Particle In Cell methods fluid simulation C++.

[![Watch the video](APIC_WAVE.PNG)](https://www.youtube.com/watch?v=wmwfUij0SW0)

## Overview:

APIC_FLUID_SIM is a fluid simulation written by me as my final major university project in C++. The fluid sim methods implemented are:

* PIC
* FLIP-PIC
* APIC

The program is split in two parts, pic-core - the fluid sim engine written in C++17 using Eigen, and point-viewer (for lack of a better name) - particle visualiser using OpenGL and SLD2 for window management.

Compatible with Ubuntu/WSLUbuntu.

## Requirements:
* c++17
* eigen3
* SDL2
* glfw
* Bazel

Build system used is BAZEL.

## Installation:

	$ sudo apt-get install bazel-2.2.0 libglfw3-dev libsdl2-dev libeigen3-dev
	$ git clone git@github.com:LHugueniot/APIC_FLUID_SIM.git

## Building:

	$ cd APIC_FLUID_SIM/pic-core
	$ bazel build //:pic_core

## Running:

	$ cd APIC_FLUID_SIM/point-viewer
	$ bazel run //:pv

Enjoy!