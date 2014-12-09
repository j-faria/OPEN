import commands
import os

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


Export('MPI', 'f2py_exec')
if MPI: Export('mpif90_exec')

## specify the gfortran executable
## --gfortran=/path/to/gfortran
AddOption('--gfortran', dest='gfortran', type='string', nargs=1, action='store',
	      help='Absolute path to Fortran compiler')


## on some machines, the gp.so extension needs to be linked with lapack explicitely
## use --gp-with-lapack if getting errors of the type
## ImportError: /data/jfaria/OPEN/OPEN/ext/gp.so: undefined symbol: dpotrf
AddOption('--gp-with-lapack', dest='gp_lapack_link', action='store_true',
	      help='some help')


SConscript(['OPEN/SConscript'])

path = 'MultiNest/SConscript'
if os.path.isfile(path):
	SConscript([path])
else:
	print "Missing the SConscript file at %s." % path
	print "Try cloning the MultiNest submodule:\n$ git submodule init\n$ git submodule update"
	Exit(1)

