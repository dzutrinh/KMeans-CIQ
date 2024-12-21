cc=gcc

all: ciq

ciq: ciq.c
	$(cc) -O2 $< -o ciq

clean:
	rm -f ciq
	
