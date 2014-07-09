CC=gcc
CFLAGS= -I/usr/include -Wall 
LDFLAGS = -L/usr/lib 
LIBS = -lsndfile -lfftw3

DEPS = fano.h
OBJ = wspr.o fano.o tab.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

wsprcan: $(OBJ) 
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS) $(LIBS)
