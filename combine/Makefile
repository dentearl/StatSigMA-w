shell:=/bin/bash -e
export SHELLOPTS=pipefail
includes := global_combine.h tree_combine.h

.PHONY = all clean test
all: ../bin/combine ../bin/getSpeciesName

../bin/combine: combine.cpp ${includes}
	mkdir -p $(dir $@)
	g++ -O2 -o $@.tmp $<
	mv $@.tmp $@

../bin/%: %.cpp
	mkdir -p $(dir $@)
	g++ -O2 -o $@.tmp $<
	mv $@.tmp $@

clean:
	rm ../bin/combine ../bin/getSpeciesName

test:
