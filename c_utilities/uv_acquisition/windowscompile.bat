ECHO "setting up variables"
CALL C:\"Program Files (x86)"\"Microsoft Visual Studio 11.0"\VC\bin\x86_amd64\vcvarsx86_amd64.bat
ECHO "running cmake"
cmake.exe -G "NMake Makefiles" ..
ECHO "running nmake"
nmake.exe /f makefile

