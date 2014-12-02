import commands

# compiling with MPI support?
out = commands.getstatusoutput('mpirun --version')
if out[0] == 0:
	MPI = True
else:
	MPI = False
	print 'Compiling MultiNest without MPI'

## get path to mpif90 executable
if MPI:
	out = commands.getstatusoutput('which mpif90')
	if out[0] == 0:
		mpif90_exec = out[1]
		mpif90_version = commands.getoutput(mpif90_exec + ' -v').split('\n')[0]
		print 'Found %s at %s' % (mpif90_version, mpif90_exec)
	else:
		raise RuntimeError('Compiling with MPI but mpif90 does not seem to be installed. Aborting!')

## get path to f2py executable
out = commands.getstatusoutput('which f2py')
if out[0] == 0:
	f2py_exec = out[1]
	f2py_version = commands.getoutput(f2py_exec + ' -v')
	print 'Found f2py version %s at %s' % (f2py_version, f2py_exec)
else:
	raise RuntimeError('f2py does not seem to be installed. Aborting!')


Export('MPI', 'mpif90_exec', 'f2py_exec')

AddOption('--gfortran', dest='gfortran', type='string', nargs=1, action='store',
	      help='Absolute path to Fortran compiler')



SConscript(['OPEN/SConscript'])
SConscript(['MultiNest/SConscript'])
