ext:
	make -C ./OPEN/ext
	@echo "All done in ext"

#nest:
#	make -C ./OPEN/multinest
#	@echo "All done in multinest"


clean:
	make clean -C ./OPEN/ext