# This makefile works on a Sandia Linux machine.
# For other platforms, you may need to adjust the "E =", "B =",
# "L =", and "W =" lines below.

E = /opt/MATLAB/R2012a
I = $E/extern/include
B = $E/bin/glnxa64
L = -L$B -leng -lmx
W = -Wl,-rpath,$B -Wl,-rpath-link,$B
CFLAGS = -O0 -g
CC = cc

dakmatlab: dakmatlab.o
	$(CC) -o dakmatlab $W dakmatlab.o $L

dakmatlab.o: dakmatlab.c
	$(CC) -c $(CFLAGS) -I$I dakmatlab.c
