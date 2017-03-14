# include the machine dependent configuration
ifneq ($(MAKECMDGOALS),clean)
  -include ../Makefile.conf
endif

.PHONY: all clean coverage

all:
ifeq ($(make_doc),yes)
	./extract
	./go-doxygen
	cp *png html/
endif

clean:
	rm -fr automatic html latex *~ links-to-user-manual.md plumedcheck.md coverage.info coverage-full.info coverage-libs.info coverage coverage-libs

coverage:
	lcov --gcov-tool gcov-4.9 -c -d ../src/ --output-file coverage.info
# remove system headers:
	lcov --gcov-tool gcov-4.9 -r coverage.info "/opt/*" "/usr/*" -o coverage-full.info
# remove third-party libraries
	lcov --gcov-tool gcov-4.9 -o coverage.info -r coverage-full.info "*/blas/*" "*/lapack/*" "*/molfile/*" "*/lepton/*"
# generate coverage
	genhtml -o coverage -t "plumed test coverage" coverage.info
# keep only third-party libraries
	lcov --gcov-tool gcov-4.9 -o coverage-libs.info -e coverage-full.info "*/blas/*" "*/lapack/*" "*/molfile/*" "*/lepton/*"
# generate coverage
	genhtml -o coverage-libs -t "plumed test coverage (other modules)" coverage-libs.info
