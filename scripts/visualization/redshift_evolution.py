import numpy as np
import h5py
import matplotlib.pyplot as plt
from parameters import Parameters
from quantities_dict import quantities_dict

plt.rc( ('axes'),labelsize=8.5)
plt.rc( ('xtick','ytick'),labelsize=8)
plt.rc("text",usetex=True)
plt.rc(("legend"),fontsize=8)

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

all_z = False
plot_type = 'log'
quantity = 'Bp'

# Opens Magnetizer output file
f = h5py.File('/data/nlfsr/magnetizer-runs/GON2_25e4_1e5_out.hdf5','r')
fi = h5py.File('/data/nlfsr/magnetizer-runs/GON2_25e4_1e5.hdf5','r')
params = Parameters(f)

rmax_rdisc = params.grid['P_RMAX_OVER_RDISK']

data = f['Output'][quantity]

ngals, ngrid, nzs = data.shape

# Target radius: half-mass radius
i_target = ngrid/rmax_rdisc
# Target radius: galaxy centre
#i_target = 1

# At the moment, no weights are being used... TODO
weights = np.ones(ngals)

results = np.empty((nzs, 3))*np.NaN # Allows invalid values to be skipped

M_min = 10**8.75
M_max = 10**9.5

M_bins = 10**np.array([8.,8.75,9.5,10.25,11.])

plt.figure(figsize=(6.5, 4.8), dpi=300) # Square

for i, (M_min, M_max) in enumerate(zip(M_bins[:-1],M_bins[1:])):
    #plt.figure()
    for iz in range(nzs-1,-1, -1):
        print 'Working on z', fi['Input']['z'][iz]

        if (iz == nzs-1) or all_z:
            # Loads to RAM the relevant part(s) of the HDF5 file (slow!)
            Mstars = fi['Input']['Mstars_disk'][:,iz] + fi['Input']['Mstars_bulge'][:,iz]

            select_mass  = Mstars > M_min
            select_mass *= Mstars < M_max

            if len(Mstars[select_mass]) == 0:
                continue

        data_selected = data[select_mass,i_target,iz]

        # Removes empty entries
        ok = data_selected>-1000

        data_selected = data_selected[ok]
        if len(data_selected ) == 0:
            continue

        # Computes percentiles
        results[iz, :] = weighted_percentile(abs(data_selected),weights[select_mass][ok])

    plt.subplot(2,2,i+1) # grid

    zs= fi['Input']['z'][:]

    plt.xlabel('z')

    if quantity in quantities_dict:
        name, units = quantities_dict[quantity]
    else:
        name, units = quantity, None

    if plot_type == 'log':
        if units:
            plt.ylabel(r'$\log({0}\,/\,{1})$'.format(name, units))
        else:
            plt.ylabel(r'$\log({0})$'.format(name))
    else:
        if units:
            plt.ylabel(r'${0}\,[{1}]$'.format(name, units))
        else:
            plt.ylabel(r'${0}$'.format(name))

    fifteen = np.log10(results[:,0])
    eightyfive = np.log10(results[:,2])
    median = np.log10(results[:,1])

    plt.plot(zs, median, marker='.', color='#d95f0e')
    plt.plot(zs, fifteen, color='#d95f0e', linestyle=':')
    plt.plot(zs, eightyfive, color='#d95f0e', linestyle=':')
    plt.fill_between(zs, fifteen, eightyfive, color='#d95f0e', alpha=0.1)
    plt.annotate(r' $ 10^{%.2f} < M_\star/{\rm M}_\odot\leq\,  10^{%.2f}$'%(
        np.log10(M_min),np.log10(M_max) ), (0.1,1.5), fontsize=9)
    plt.xlim([0,4.9])
    plt.ylim([-4.5,2.5])

plt.subplots_adjust(left=0.07,
                    right=0.99,
                    bottom=0.075,
                    top=0.98)

if all_z:
  plt.savefig('/home/nlfsr/zevol_halfmass.pdf')
else:
  plt.savefig('/home/nlfsr/zevol_halfmass_z0.pdf')
