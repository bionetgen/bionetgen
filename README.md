# Bionetgen

Bionetgen is a brief C++ tool to create fiber networks that follow the fiber length, valency and cosine distribution of actual collagen gels from confocal microscope images. Our algorithm closely follows the approach of [1] using the stochastic optimization method of simulated annealing.

[1] S.B. Lindstrom, D.A. Vader, A. Kulachenko, D.A. Weitz, Biopolymer network geometries: Characterization,
regeneration, and elastic properties. Phys. Rev. E - Stat. Nonlinear, Soft Matter Phys. 82(5), 2 (2010).


## Usage

    $ ./build/voronoi config.json


## Configuration

Use the config.json file to set parameters for the algorithm.

The *main* operation mode can be selected by changing the values of *generate* and *simulate* in the config file. If only generate is true, a voronoi geometry will be generated and saved as the output. If simulate is also true, simulated annealing will additionally be performed before the geometry is saved to the output files. If only simulate is true, a previously saved geometry will be read from the input files, simulated annealing will be perfomed and the result will be saved to the output files.


## Building

### Building on Linux

Assumes you have `cmake`, `make`, and a C/C++ compiler installed.

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make
    $ cd ..


### Building on Windows

These build instructions use the MinGW64 part of [MSYS2](https://www.msys2.org/), which provides the necessary tools (i.e. `cmake`, `make`, `g++`, and `boost`). It generates a standalone executable that only depends on Windows' core libraries.

The reason this guide uses MinGW64 under MSYS2, despite MSYS2 generally being considered to be more modern, is because MinGW can be used to create standalone executables that ultimately do not depend on either MSYS2 or MinGW64 to run them.

- Install [MSYS2](https://www.msys2.org/)
- Clone this repository onto your computer
- Boot into an `MSYS2_MinGW64` terminal
- Go through these commands:

```bash
# install `cmake`, `make`, `gcc`, `g++`, and `boost`
pacman -S mingw-w64-x86_64-cmake make mingw-w64-x86_64-gcc mingw-w64-x86_64-boost

# change to the cloned repository
cd /C/Users/user/Desktop/bionetgen

# build with static linking
mkdir build
cd build
cmake .. -DCMAKE_EXE_LINKER_FLAGS="-static"
cmake --build .

# check voronoi.exe exists, and that it only depends on Windows libraries
#
# (i.e. every library it lists should be in /C/Windows)
ldd voronoi.exe
```

#### Developing on Windows (Visual Studio 2022)

- Confirm you can build `voronoi.exe` manually using MSYS2 etc. (above)
- Install a debugger (`gdb`) in your MSYS2 terminal:

```bash
pacman -S mingw-w64-x86_64-gdb
```

- In Visual Studio 2022, open `bionetgen` as a folder project ("Open as Local Folder" from the splash screen)
- Visual Studio should detect that it's a CMake project and suggest `Open CMake Settings Editor`
- Open the CMake settings editor and save `Ctrl+S`, which should create a `CMakeSettings.json` file
- Right-click the `CMakeSettings.json` file and `Open With` the JSON editor (rather than the default GUI editor)
- Replace the original content with this custom JSON (it configures Visual Studio to use MSYS2 to build the C++)

```json
{
  "configurations": [
    {
      "environments": [
        {
          "MINGW64_ROOT": "C:/msys64/mingw64",
          "BIN_ROOT": "${env.MINGW64_ROOT}/bin",
          "FLAVOR": "x86_64-w64-mingw32",
          "TOOLSET_VERSION": "12.2.0",
          "INCLUDE": "${env.MINGW64_ROOT}/include/c++/${env.TOOLSET_VERSION};${env.MINGW64_ROOT}/include/c++/${env.TOOLSET_VERSION}/${env.FLAVOR};${env.MINGW64_ROOT}/include/c++/${env.TOOLSET_VERSION}/backward;${env.MINGW64_ROOT}/lib/gcc/${env.FLAVOR}/${env.TOOLSET_VERSION}/include;${env.MINGW64_ROOT}/lib/gcc/${env.FLAVOR}/${env.TOOLSET_VERSION}/include-fixed;${env.MINGW64_ROOT}/include",
          "MINGW_PREFIX": "C:/msys64/mingw64",
          "MINGW_CHOST": "x86_64-w64-mingw32",
          "MINGW_PACKAGE_PREFIX": "mingw-w64-x86_64",
          "MSYSTEM": "MINGW64",
          "MSYSTEM_CARCH": "x64_64",
          "MSYSTEM_PREFIX": "${env.MINGW64_ROOT}",
          "MSYSTEM_CHOST": "x86_64-w64-mingw32",
          "SHELL": "${env.MINGW_PREFIX}/../usr/bin/bash",
          "TEMP": "${env.MINGW_PREFIX}/../tmp",
          "TMP": "${env.TEMP}",
          "PATH": "${env.MINGW_PREFIX}/bin;${env.MINGW_PREFIX}/../usr/local/bin;${env.MINGW_PREFIX}/../usr/bin;${env.MINGW_PREFIX}/../bin;${env.PATH}",
          "environment": "mingw_64_custom"
        }
      ],
      "name": "Mingw64-Debug",
      "generator": "Ninja",
      "configurationType": "Debug",
      "inheritEnvironments": [ "mingw_64_custom" ],
      "cmakeCommandArgs": "",
      "buildCommandArgs": "-v",
      "ctestCommandArgs": "",
      "variables": [
        {
          "name": "CMAKE_C_COMPILER",
          "value": "${env.BIN_ROOT}/gcc.exe",
          "type": "STRING"
        },
        {
          "name": "CMAKE_CXX_COMPILER",
          "value": "${env.BIN_ROOT}/g++.exe",
          "type": "STRING"
        }
      ]
    }
  ]
}
```

- Save the JSON, which should cause Visual Studio to reconfigure the project
- Reconfiguration should succeed. The `Output` tab (bottom of Visual Studio)
  will print something like "Configuration Complete". If it doesn't, there's
  probably something wrong in the JSON. Good guesses are:

    - `MINGW64_ROOT` is wrong, because you installed MSYS2 somewhere else. Find
       where it's installed and change it in the JSON to match the actual location
    - `TOOLSET_VERSION` is wrong, because your MSYS2 installed a slightly different
      version of C++. This makes later strings like (e.g.) `${env.MINGW64_ROOT}/include/c++/${env.TOOLSET_VERSION}`
      incorrect. The solution is to browse through your MSYS2 install to find what
      version string to use (look at how it's used in the example)

- In Visual Studio, at the top of the UI where it says "Select Startup Item", click
  the little dropdown arrow and select `voronoi.exe`. This will cause (e.g.) `Ctrl+B`
  to build `voronoi.exe` and (e.g.) `F5` to build+debug `voronoi.exe`

- If you want to run `voronoi.exe` with arguments (e.g. a file, input configuration, etc.)
  then you will need to tell Visual Studio which arguments to use.

  - Stackoverflow guide: https://stackoverflow.com/questions/30104520/adding-command-line-arguments-to-project
  - Microsoft guide: https://devblogs.microsoft.com/cppblog/using-mingw-and-cygwin-with-visual-cpp-and-open-folder/
  - Build environments: https://learn.microsoft.com/en-us/cpp/build/cmake-predefined-configuration-reference?view=msvc-170

- Effectively, you need to right-click `CMakeLists.txt` then `Add Debug Configuration`. You might
  get lucky and find that Visual Studio suggests "Run with debugger mingw64", which will generate
  a mostly-correct configuration. Otherwise, you may need to `Open Debug and Launch Settings`
  to open the JSON followed by writing the appropriate JSON. Here is what worked for me. It generates
  a `voronoi.exe` startup target *and* a `voronoi.exe (gdb)` target that attaches a debugger to the
  program when it boots:

```json
{
  "version": "0.2.1",
  "defaults": {},
  "configurations": [
    {
      "type": "default",
      "project": "CMakeLists.txt",
      "projectTarget": "voronoi.exe",
      "name": "voronoi.exe",
      "args": [ "C:\\some\\path\\to\\config.json" ]
    },
    {
      "type": "cppdbg",
      "project": "CMakeLists.txt",
      "projectTarget": "voronoi.exe",
      "cwd": "${workspaceRoot}",
      "program": "${debugInfo.target}",
      "name": "voronoi.exe (gdb)",
      "args": [ "C:\\some\\path\\to\\config.json" ],
      "MIMode": "gdb",
      "miDebuggerPath": "C:\\msys64\\mingw64\\bin\\gdb.exe",
      "externalConsole": true
    }
  ]
}
```

### Publications

Please cite the following paper for acknowledging this software contribution:

Eichinger JF, Grill MJ, Davoodi Kermani I, Aydin RC, Wall WA, Humphrey JD, Cyron CJ, A computational framework for modeling cell-matrix interactions in soft biological tissues. *Biomechanics and Modeling in Mechanobiology*. 2021. doi.org/10.1007/s10237-021-01480-2

### Parameters

*seed*: The seed for the random number generator (can be any integer number)

*particles*: The number of Voronoi particles (can be any number greater 5)

*input-prefix*: The path prefix of the input files if simulated annealing is done on a previously generated geometry (relative to the working directory, optional)

*output-prefix*: The path prefix for the output files (relative to the working directory, not the )

*generate*: If a voronoi geometry should be generated (can be true or false)

*simulate*: If simulated annealing should be performed on a voronoi geometry (can be true or false)

*box-size*: The size of the box containing the voronoi geometry to be generated (must be an array of 3 positive values)

*box-origin*: The origin of the box (must be an array of 3 values)

### *simulated-annealing*: Options for simulated annealing (if needed)

*mode*: Simulated annealing mode (can be 1, 2 or "both")

*max-iter*: Maximum number of iterations  (can be any positive integer)

*max-subiter*: Maximum number of sub-iterations (can be any positive integer)

*weight-line*: Line weight (can be any positive value)

*weight-cosine*: Cosine weight (can be any positive value)

*tolerance*: Tolerance (can be any positive value)

*temperature-initial*: Initial temperature (can be any positive value)

*temperature-decay-rate*: Temperature decay rate  (can be any positive value)

*max-movement-frac*: Maximum movement fraction (can be any value between 0 and 1)

*screen-output-every*: After how many iterations screen output is produced (can be any positive integer)

*num-bins-length*: Number of bins per length (can be any positive integer)

*num-bins-cosine*: Number of bins per cosine (can be any positive integer)

### License

*bionetgen* is published under the BSD 3-Clause License.

### Acknowledgements

This project uses the library voro++ by Chris Rycroft from University of California, through Lawrence Berkeley National Laboratory, for the generation of the voronoi geometry, which can be downloaded from http://math.lbl.gov/voro%2B%2B/.

