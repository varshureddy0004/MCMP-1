CC_GNU = gcc
CC_INTEL = icc
CFLAGS = -O3 -fopenmp
INTEL_FLAGS = -O3 -qopenmp

all: ssi mmi ssomp mmomp issomp immomp

ssi: ssInsertion.c
	$(CC_GNU) $(CFLAGS) ssInsertion.c -o ssi.exe -lm

mmi: mmInsertion.c
	$(CC_GNU) $(CFLAGS) mmInsertion.c -o mmi.exe -lm

ssomp: ompssInsertion.c
	$(CC_GNU) $(CFLAGS) ompssInsertion.c -o ssomp.exe -lm

mmomp: ompmmInsertion.c
	$(CC_GNU) $(CFLAGS) ompmmInsertion.c -o mmomp.exe -lm

issomp: ompcInsertion.c
	$(CC_INTEL) $(INTEL_FLAGS) ompcInsertion.c -o issomp.exe

immomp: ompfInsertion.c
	$(CC_INTEL) $(INTEL_FLAGS) ompfInsertion.c -o immomp.exe

clean:
	rm -f *.exe

.PHONY: all clean