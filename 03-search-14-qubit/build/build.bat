@echo off
setlocal
call "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cl /nologo /std:c++20 /O2 /Oi /Ot /GL /EHsc /W4 /permissive- /DNDEBUG /arch:AVX2 stab14.cpp /Fe:stab14.exe /link /LTCG
