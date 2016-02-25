""" Contains routines which allow extracting galform data from an output hdf5
    file and repack it in a more usable way. """
import numpy as N
import h5py
import os
from numpy.random import random, shuffle
#from empirical_disks import empirical_disk
from data_dict_util import filter_dictionary_inplace, filter_dict

def get_parameter_values(model_dir):
    """read the used_paramaters file in the directory and list the
    parameter values.  Note that values are returned as strings"""

    if os.path.isfile( model_dir+"/used-parameters.bz2" ):
        os.system("bunzip2 -f "+model_dir+"/used-parameters.bz2")
    #(params, values) = N.genfromtxt( model_dir + "/used-parameters", \
    #                          unpack=True, usecols=(0,1), delimiter="=",\
    #                          dtype=None, autostrip=True)

    # it is more useful to create a dictionary with the parameter as the key
    d = dict( N.genfromtxt( model_dir + "/used-parameters", \
                              unpack=False, usecols=(0,1), delimiter="=",\
                              dtype=None, autostrip=True))
    for p in d.keys():            # where possible convert strings to values.
        try:
            d[p] = N.float64( d[p] )
        except:
            pass

    return d

def read_time_data(model_dir, maximum_final_B_over_T=0.5, return_data_dict=True,
                   minimum_final_stellar_mass=1e6, minimum_final_gas_mass=1e5,
                   number_of_galaxies=100, empirical_disks=True, ivol_dir='',
                   minimum_final_disk_size=0.0,datasets=None):
    """Reads data from the galaxies.hdf5 file inside model_dir/ivol_dir (for
       Monte-Carlo runs, ivol_dir can be left blanck). The option
       number_of_galaxies sets the _approximate__ number of galaxies in the
       output (which are drawn, randomly, from the sample, taking into account
       halo press-schechter weights. Other options allow to select by minimum
       mass of the final galaxy.

       Returns a dictionary of dictionaries containing times as keys in the
       first level, dataset names in the second level and the arrays
       associated with each dataset in the last level.
       Redshift and time information can be obtained from the special keys: zout
       and tout.  """


    if datasets == None:
        # Which data sets will be read from the output
        datasets = (
                    # Essential physical properties
                    'vdisk', 'mstars_disk', 'mstars_bulge', 'mcold', 'mhot',
                    'mcold_burst','rdisk', 'rbulge', 'mstardot', 'mhalo',
                    'is_central', 'mstardot_average',
                    'vbulge','vchalo','halo_r_virial','strc',
                    # Old-style indexing
                    'ident_final', 'ident_next_output', 'index',
                    # New style indexing and extras
                    'FirstProgenitorID',
                    'EndMainBranchID','DescendantID',
                    'GalaxyID', 'LastProgenitorID',
                   )

    # Gets the parameters (later this may be used...)
    params = get_parameter_values(model_dir+'/'+ivol_dir)
    h0     = params['h0']

    # Opens the hdf5 file
    if os.path.isfile( model_dir+'/'+ivol_dir+"/galaxies.hdf5.bz2" ):
        os.system("bunzip2 -f "+model_dir+'/'+ivol_dir+"/galaxies.hdf5.bz2")
    print model_dir+'/'+ivol_dir+"/galaxies.hdf5"
    f = h5py.File(model_dir+'/'+ivol_dir+"/galaxies.hdf5")

    # Gets the redshifts and times
    zout_array = f["Output_Times/zout"][:]
    tout_array = f["Output_Times/tout"][:]
    z_indices = N.arange(len(zout_array))

    data_dict = {}

    for i, zidx in enumerate(z_indices[:]):
        t = tout_array[zidx]
        print('Working on redshift {0}, time {1}'.format(zout_array[zidx],t))

        data_dict[t] = {}
        output_id = "Output%.3i" % (zidx+1)   # outputs labeled Fortran style.

        # read the data from the hdf5 file.
        for k in datasets:
            data_dict[t][k] = f[output_id][k]

        if "galaxy_weight" not in f[output_id]:
            print('Adding weights, jm and jtree (this may take some time).')
            # Weights are associated with haloes.
            # Need to assign a weight to each galaxy.
            # The same goes with jm and jtree
            halo_weight = f[output_id+"/Trees/weight"]
            halo_ngals = f[output_id+"/Trees/ngals"]
            halo_jm = f[output_id+"/Trees/jm"]
            halo_jtree = f[output_id+"/Trees/jtree"]

            if ((N.sum(halo_ngals) != len(data_dict[t]['mstars_disk']) or 
                 len(halo_weight) != len(halo_ngals))):
                print "ERROR: not enough galaxies in the file"
                exit()
            
            # count galaxies, assiging weight from their halo weight.
            # initial weights are -1
            galaxy_weight = N.zeros( len(f[output_id+"/mstars_disk"][:]) )-1.
            jtree = N.zeros( len(f[output_id+"/mstars_disk"][:]) ) -1.
            jm = N.zeros( len(f[output_id+"/mstars_disk"][:]) ) -1.

            igal = 0
            for ihalo in range( len(halo_ngals) ):
                if (halo_ngals[ihalo] > 0):
                    galaxy_weight[igal:igal+halo_ngals[ihalo]] = halo_weight[ihalo]
                    jtree[igal:igal+halo_ngals[ihalo]] = halo_jtree[ihalo]
                    jm[igal:igal+halo_ngals[ihalo]] = halo_jm[ihalo]
                    igal = igal+halo_ngals[ihalo]

            f[output_id+"/galaxy_weight"] = galaxy_weight
            f[output_id+"/galaxy_jtree"] = jtree
            f[output_id+"/galaxy_jm"] = jm
            f.flush()
        else:
            galaxy_weight = f[output_id+"/galaxy_weight"]

        data_dict[t]['weight'] = galaxy_weight
        
        print 'Galaxies read:',len(galaxy_weight)
        if i==0:
            print "Pre-sampling {0} galaxies based on weight".format(
                                                        10.0*number_of_galaxies)

            randm = random( len(galaxy_weight) )
            # select a number of galaxies based on weight
            if 10*number_of_galaxies < len(galaxy_weight) :
                rfac = (10.0*number_of_galaxies) / N.sum(galaxy_weight)
                ok = ( randm / data_dict[t]['weight'] < rfac )

                filter_dictionary_inplace(ok,data_dict[t])

            data_dict[t]['BoT'] = (data_dict[t]['mstars_bulge'][:]+
                                    data_dict[t]['mcold_burst'][:]) / (
                                        data_dict[t]['mstars_disk'][:]+
                                        data_dict[t]['mcold'][:]+
                                        data_dict[t]['mstars_bulge'][:]+
                                        data_dict[t]['mcold_burst'][:] )

            # Selects galaxies with a minimum stellar mass, gass mass and
            # bulge to total mass ratio
            ok  = data_dict[t]['BoT'][:] < maximum_final_B_over_T
            ok *= data_dict[t]['mstars_disk'][:] > minimum_final_stellar_mass
            ok *= data_dict[t]['mcold'][:] > minimum_final_gas_mass
            ok *= data_dict[t]['rdisk'][:] > minimum_final_disk_size
            filter_dictionary_inplace(ok,data_dict[t])

            print('Number of galaxies after filtering: {0}'.format(
                                                   len(data_dict[t]['weight'])))
            print("Sampling {0}  galaxies based on weight".format(
                                                            number_of_galaxies))
            randm = random( len(data_dict[t]['weight']) )
            # selects a number of galaxies based on weight
            rfac = number_of_galaxies / N.sum(data_dict[t]['weight']) 
            ok = ( randm / data_dict[t]['weight'] < rfac )
            
            print('Actual number of selected galaxies: {0}'.format(len(ok[ok])))
            filter_dictionary_inplace(ok,data_dict[t])

            print('Number of galaxies after sampling: {0}'.format(
                                                   len(data_dict[t]['weight'])))

            # Stores indexing information and reference ID
            GalaxyID_target    = data_dict[t]['FirstProgenitorID']
            data_dict[t]['ID']= data_dict[t]['GalaxyID']
            previous_ID = data_dict[t]['ID']
            previous_t = t
        else: 
            # selects only the most massive ascendents of the gals in
            # the previous t
            filt = N.zeros(data_dict[t]['GalaxyID'].shape, dtype=bool)
            data_dict[t]['ID'] = N.zeros(data_dict[t]['GalaxyID'].shape) * N.NaN

            for target_gid, ID in zip(GalaxyID_target,previous_ID):
                # Finds the index of a target galaxy
                match = target_gid == data_dict[t]['GalaxyID']

                if match.any():
                    # If there is a match, updates filter
                    filt +=  match
                    # Includes the z=0 ID in the dictionary
                    data_dict[t]['ID'][match] = ID

            # Applies the filter
            filter_dictionary_inplace(filt,data_dict[t])

            print('Number of galaxies after sampling: {0}'.format(
                                                      len(data_dict[t]['ID'])))

            # Creates a list of missing galaxies

            GalaxyID_target = data_dict[t]['FirstProgenitorID'][:]
            previous_ID = data_dict[t]['ID']
            previous_t = t



    if empirical_disks: # TODO not working yet
        pass
        # Deep copy of data_dict_z0, to be used in the empirical_disk function
        #data_dict_z0 = data_dict[first_t].copy()
        #for t in data_dict:
            #data_dict[t] = empirical_disk(data_dict[t],data_dict_z0)

    data_dict['tout'] = tout_array
    data_dict['zout'] = zout_array

    return data_dict

def plot_mass_evolution(model_dir, gtype='all'):
    """ Tests the module ploting the time evolution of galaxy mass.
        gtype can be either: all, central/cen or sat/satellite."""
    import pylab as P

    # Will color them by mass
    cmap = P.cm.get_cmap('YlGnBu',200)

    data_dict = read_time_data(model_dir,
                                maximum_final_B_over_T=1.0,
                                minimum_final_stellar_mass=1e10,
                                minimum_final_gas_mass=1e7,
                                number_of_galaxies=100,
                                empirical_disks=False,
                                ivol_dir='ivol0')

    ts = sorted(data_dict['tout'])[::-1]
    IDs = data_dict[ts[0]]['ID']
    if gtype == 'all':
        pass
    elif 'sat' in gtype:
        IDs = IDs[data_dict[ts[0]]['is_central']==0]
    elif 'cen' in gtype:
        IDs = IDs[data_dict[ts[0]]['is_central']==1]
    else:
        exit('plot_mass_evolution: Error! Unrecongnized gtype.')

    shuffle(IDs)
    IDs = IDs[:10]

    gals_dict = dict()
    max_mass=None

    for ID in sorted(IDs):
        gals_dict[ID] = []
        for t in ts:
            select = data_dict[t]['ID'] == ID
            mass  = (data_dict[t]['mcold'][select] +
                    data_dict[t]['mstars_disk'][select] +
                    #data_dict[t]['mhot'][select] +
                    data_dict[t]['mstars_bulge'][select]+
                    data_dict[t]['mcold_burst'][select])

            if len(mass)==1:
                if max_mass==None:
                    max_mass = mass[0]
                    min_mass = mass[0]
                max_mass = max(mass[0], max_mass)
                min_mass = min(mass[0], min_mass)
                gals_dict[ID].append(mass[0])
            else:
                gals_dict[ID].append(N.NaN)

    for ID in IDs:
        c=cmap( (gals_dict[ID][0]-min_mass)/(max_mass-min_mass) )
        P.plot(ts,N.array(gals_dict[ID])/1e10, marker='.', color=c)
    P.xlabel(r'$t\,[{{\rm Gyr}}]$')
    P.ylabel(r'$M\,/\, (10^{10}\,{{\rm M}}_\odot)$')
    P.show()

if __name__ == "__main__"  :
    model_dir = 'test_SAM_output'
    plot_mass_evolution(model_dir, gtype='cen')
