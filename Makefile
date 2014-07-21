NO_COLOR=\033[0m
OK_COLOR=\033[92m
ERROR_COLOR=\033[91m
WARN_COLOR=\033[93m

OK_STRING=$(OK_COLOR)[OK]$(NO_COLOR)
ERROR_STRING=$(ERROR_COLOR)[ERRORS]$(NO_COLOR)
WARN_STRING=$(WARN_COLOR)[WARNINGS]$(NO_COLOR)

all: ext nest

ext: 
	@echo "Compiling Fortran extensions..."
	@make -C ./OPEN/ext >> compile.out
	@echo "OPEN extensions  --  $(OK_STRING)"

.PHONY: MultiNest
MultiNest: 
	@echo "Compiling MultiNest..."
	@make -C ./MultiNest/MultiNest_v3.7 --quiet 
	@make libnest3.so -C ./MultiNest/MultiNest_v3.7 --quiet 
	@echo "MultiNest  --  $(OK_STRING)"

nest: MultiNest
	@echo "Compiling OPEN-MultiNest interface..."
	@make clean -C ./OPEN/multinest --quiet 
	@make fortran -C ./OPEN/multinest --quiet 
	@echo "OPEN <-> MultiNest  --  $(OK_STRING)"


clean:
	make clean -C ./OPEN/ext
	make clean -C ./OPEN/multinest
	make clean -C ./MultiNest/MultiNest_v3.7
