# EFMCUtils
These are soft-coded macros that I generally use in my ROOT sesions.

To use the macros in your sesions follow these steps in you terminal:

1. Go to the folder containing EFMCUtils.C . This folder will be called "source"
2. Open a ROOT session
3. Type:
```
.L EFMCUtils.C++
```
This will compile the macros into a library that ROOT can read
4. ROOT checks several options and folder when it is opened. Therefore, we have to specify to it where is located this new library containing our macros. For this matter we modify (or create) the .rootrc file in the ${HOME} folder. Add to this file:
```
Unix.*.Root.DynamicPath:     .:${ROOTSYS}/lib:${HOME}/EFMCUtils/lib:
Unix.*.Root.MacroPath:       .:${HOME}/EFMCUtils:

Root.ShowPath: false
```
In the folder containing EFMCUtils.C add a rootlogon.C file with the following line:
```
gROOT->LoadMacro("EFMCUtils.C");
```
This should be enough. In case ROOT do not recognizes the macros, add to the file .rootlogon.C, which is usually located at ${HOME}, the following line:
```
gROOT->ProcessLine(".L ${HOME}/EFMCUtils/EFMCUtils.C+");
```
