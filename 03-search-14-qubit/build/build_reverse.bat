@echo off
setlocal
call "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cl /nologo /std:c++20 /O2 /Oi /Ot /GL /EHsc /W4 /permissive- /D_CRT_SECURE_NO_WARNINGS /DNDEBUG /arch:AVX2 reverse_8w10_3w8.cpp /Fe:reverse_8w10_3w8.exe /link /LTCG
