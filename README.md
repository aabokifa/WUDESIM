# WUDESIM
A special water quality simulation module for the dead-end branches of drinking water distribution systems. The module is an add-on to the existing EPANET water quality simulation module.

To build using CMake and Visual Studio on Windows: 

- Make sure that CMake is installed and accessible from the system PATH. You can check by typing "cmake" at the command prompt which should result in a message from the program, and not an error like " 'cmake' is not recognized ".

- Open command prompt and navigate to the project root directory. You can easily do this in Windows by holding down the Shift key while right clicking in the directory folder.

Type in the following commands:

> mkdir build

> cd build

> cmake ..

CMake will identify the installed compilers on your machine, and automatically build the project for you using the latest version of Visual Studio. 

- In the build directory, a Visual Studio solution file named "WUDESIMmake.sln" is created. To build the executable file "WUDESIM.exe", open the solution file in Visual Studio, right click on WUDESIM in the solution explorer and choose Build. The executable file is created in the /build/src/Debug directory.

- Alternatively, you can directly generate the executable by typing the following command after building the project with CMake:
>msbuild /p:Configuration=Release ALL_BUILD.vcxproj

The executable will be generated in the /build/src/Release directory.

- To run the executable from the command line:
>WUDESIM.exe EPANET.inp EPANET.rpt WUDESIM.inp

- In the /Sample directory you will find three sample files (Net2.inp, Net2.rpt, and WUDESIM.inp)to run the program with.
