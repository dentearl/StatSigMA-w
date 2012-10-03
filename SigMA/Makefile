shell:=/bin/bash -e
export SHELLOPTS=pipefail
includes := alpha.h blast.h block.h global.h protmat.h pvalue.h tree.h common.h

CPPFLAGS = -O3 -g -Wall -Werror -pedantic

.PHONY =  all clean test
all: ../bin/SigMA

../bin/SigMA: SigMA.cpp ${includes}
	mkdir -p $(dir $@)
	g++ ${CPPFLAGS} -o $@.tmp $<
	mv $@.tmp $@

clean:
	rm ../bin/SigMA

test:
