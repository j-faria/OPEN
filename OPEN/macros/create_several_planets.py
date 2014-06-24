# 'run -i' this script
magic = get_ipython().magic

P = 22
e = 0.5
k = 43
N = 123
f = 'testedata%d.dat' % 2


param = 'p(%f) e(%f) k(%f) N(%d) file(%s)' % (P, e, k, N, f)
print param

magic('create ' + param)