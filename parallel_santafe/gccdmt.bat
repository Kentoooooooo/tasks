@echo off

REM gcc -O3 -c -DDSFMT_MEXP=19937 -o dSFMT.o dSFMT.c
REM ar rcs dSFMT.a dSFMT.o

SET LIB=C:/Users/18_LaserNetwork/Documents/Num-Calc/rc/multi_fb/lib
SET HED=C:/Users/18_LaserNetwork/Documents/Num-Calc/rc/multi_fb/inc

REM gcc -Wall -O2 -DDSFMT_MEXP=19937 -I%HED% -L%LIB% %1 -ldSFMT
gcc -O2 -DDSFMT_MEXP=19937 -I%HED% -L%LIB% %1 -ldSFMT
