import numpy as np
import healpy as hp

nside = 10  #nDirections = 12 * nside^2
center = [0.48030293, 0.49758723, 0.5158643]
radius = 0.00358
nBins = 300
lmin = 1299
lmax = 1307.5
nvec = hp.nside2npix(nside)

m = np.arange(nvec)
k = [[0]*3]*nvec
for i in range (0,nvec) :
  k[i] = hp.pix2vec(nside,i)


f = open('params_mock_1200.dat','w')

for i in range(nvec):
    f.write('%07.5f %07.5f %07.5f \n'%(k[i][0],k[i][1],k[i][2]))
    f.write('%07.5f, %07.5f, %07.5f \n'%(center[0],center[1],center[2]))
    f.write('0.0 \n')
    f.write('%i %07.5f %6.1f %6.1f \n'%(nBins, radius, lmin, lmax))
    f.write('0 0.0 \n')
    f.write('0 0 0 0 0. \n')
    f.write('\n')



