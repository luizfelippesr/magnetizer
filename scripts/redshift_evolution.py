import numpy as np
import h5py
import matplotlib.pyplot as plt

def weighted_percentile(data, weights, percentiles=[15,50,85]):
      """
      Computes percentiles
      Input: data -> 1D-numpy array containing the data
             weights -> 1D-numpy array containing the weights
             percentiles -> a sequence of the percentiles to be returned
                            Default: [15,50,85]
      Output: an array containing the requested percentiles

      """
      # Sorts the arrays
      ind = np.argsort(data)
      sorted_data = data[ind]
      sorted_weights = weights[ind]

      # Computes all the percentiles
      P = 100.0*(sorted_weights.cumsum()-0.5*sorted_weights)/sorted_weights.sum()

      # Interpolates the requested percentiles
      return np.interp(percentiles, P, sorted_data)

# Opens Magnetizer output file
f = h5py.File('/home/nlfsr/magnetizer_runs/GON_test_1500_output.hdf5','r')

# TODO This should be set using the data in the HDF5 file
rmax_rdisc = 2.25

Bphi = f['Output']['Bp']

ngals, ngrid, nzs = Bphi.shape

# Target radius: half-mass radius
i_target = ngrid/rmax_rdisc
# Target radius: galaxy centre
#i_target = 1

# At the moment, no weights are being used... TODO
weights = np.ones(ngals)

results = np.empty((nzs, 3))
for iz in range(nzs):
    # Loads to RAM the relevant part of the HDF5 file (slow!)
    Bphi_selected = Bphi[:,i_target,iz]
    # Removes empty entries
    ok = Bphi_selected>-1000
    # Computes percentiles
    results[iz, :] = weighted_percentile(Bphi_selected[ok],weights[ok])

zs= f['Input']['z'][:]

plt.plot(zs, results[:,1], marker='.')
plt.plot(zs, results[:,0])
plt.plot(zs, results[:,2])
plt.show()
