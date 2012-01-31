shell:=/bin/bash -e
export SHELLOPTS=pipefail
includes := global_combine.h tree_combine.h

CPPFLAGS = -O2 -g -Wall -Werror

.PHONY = all clean test
all: ../bin/combine ../bin/getSpeciesName

../bin/%: %.cpp ${includes}
	mkdir -p $(dir $@)
	g++ ${CPPFLAGS} -o $@.tmp $<
	mv $@.tmp $@

clean:
	rm ../bin/combine ../bin/getSpeciesName

test:
