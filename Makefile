CC      ?= gcc
CFLAGS  ?= -O2 -Wall -Wextra -std=c99
PY      ?= python

.PHONY: all test-c reproduce milp clean

all: test-c reproduce

## Build the C implementation and check it against the official NSA test vector.
test-c:
	$(CC) $(CFLAGS) impl/aradi.c impl/test_vector.c -o impl/test_vector
	./impl/test_vector

## Headline experiment: cube attack on modified 6-round ARADI, 100 random keys.
reproduce:
	cd python && $(PY) run_multicube_100.py

## MILP verification of the AABB cube-distinguisher.
milp:
	cd python && $(PY) run_6round.py

clean:
	rm -f impl/test_vector impl/test_vector.exe python/*_results*.txt
