CC=gcc
CFLAGS= -I/usr/include -Wall -O2
LDFLAGS = -L/usr/lib 
LIBS = -lsndfile -lfftw3 -lm

DEPS = fano.h
OBJ = wspr.o fano.o tab.o nhash.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

k9an-wsprd: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS) $(LIBS)
