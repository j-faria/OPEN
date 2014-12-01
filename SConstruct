import commands

# compiling with MPI support?
out = commands.getstatusoutput('mpirun --version')
if out[0] == 0:
	MPI = True
else:
	MPI = False
	print 'Compiling MultiNest without MPI'

Export('MPI')

AddOption('--gfortran', dest='gfortran', type='string', nargs=1, action='store',
	      help='Absolute path to Fortran compiler')



SConscript(['OPEN/SConscript'])
# SConscript(['MultiNest/SConscript'])
