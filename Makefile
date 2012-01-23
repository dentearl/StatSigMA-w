shell:=/bin/bash -e
export SHELLOPTS=pipefail
packages = SigMA combine

.PHONY: all clean test

all: $(foreach f,${packages},$f.all)

clean: $(foreach f,${packages},$f.clean)

test: $(foreach f,${packages},$f.test)

%.all:
	cd $* && make all

%.clean:
	cd $* && make clean

%.test:
	cd $* && make test
