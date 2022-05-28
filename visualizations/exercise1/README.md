# Computer Graphics 2

This repository contains a C++ skeleton for the Computer Grapics 2 exercises at TU Berlin.

## Requirements

- Development tools (compiler, dependencies, ...)
    - Ubuntu/Debian:
      ```bash
      sudo apt install build-essential
      sudo apt install xorg-dev libglu1-mesa-dev freeglut3-dev mesa-common-dev # Polyscope dependencies
      ```
    - Windows: Install [Visual Studio 2022](https://visualstudio.microsoft.com/de/thank-you-downloading-visual-studio/?sku=Community&channel=Release&version=VS2022&source=VSLandingPage&cid=2030&passive=false) with C++ desktop workload (Community edition is free)

- git
    - Ubuntu/Debian: 
      ```bash
      sudo apt install git
      ```
    - Windows: Install [Git for Windows](https://git-scm.com/download/win)


- CMake >= 3.14
    - Ubuntu/Debian: 
      ```bash
      sudo apt install cmake
      ```
    - Windows: Install [CMake from the website](https://cmake.org/download/)

## Compilation

In the root directory, run the following commands

Linux Makefiles:

```bash
cmake . -B build -DCMAKE_BUILD_TYPE=RelWithDbInfo # BUILD_TYPE can also be `Release` or `Debug`
cmake --build build --parallel
```

Visual Studio:

```bash
cmake . -B build
cmake --build build --parallel
```

The commands compile all targets by default. If you want to compile a specific target (e.g. `ex1`) use

```bash
cmake --build build --target ex1 --parallel
```

The compiled executables are located in the `build/bin` subdirectory.

