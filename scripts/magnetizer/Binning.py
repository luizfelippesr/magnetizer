import numpy as np
import astropy.units as u


class BinningObject(object):
    """
    The attributes contain the mass bins and filters which can be applied to
    select specific mass bins.
    """
    def __init__(self, magnetizer_run, z=0.0, bins=None,
                 bin_array=None):

        self.run = magnetizer_run # Reference to run object (to avoid mistakes)
        if bins is None:
            if bin_array is None:
                raise ValueError
            M_bins = bin_array
            bins = []
            for m_m, m_p in zip(M_bins[:-1],M_bins[1:]):
                bins.append((m_m, m_p))

        self.masks = [None]*len(self.bins)
        self.bins = tuple(bins)
        self._quantitytype = None

    def _compute_bin_filter(self,quantity, bin_interval):
        print quantity
        select  = quantity >  bin_interval[0]
        select *= quantity <= bin_interval[1]
        return select

    def _update_masks(self,quantity):
        for i, interval in enumerate(self.bins):
            mask = self._compute_bin_filter(quantity, interval)
            if self.masks[i] is None:
               self.masks[i] = mask
            else:
               self.masks[i] *= mask

    def __repr__(self):
        return_str = '[BinningObject]\n'
        return return_str + self.__str__()

    def __str__(self):
        return_str = '{} bins:\n'.format(self._quantitytype)
        for (b1, b2) in self.bins:
            return_str += '{0}   {1}\n'.format(b1, b2)
        return return_str


class MassBinningObject(BinningObject):
    def __init__(self, magnetizer_run, z=0.0, bins=None, stellar_mass=True,
                 gas_mass=False, include_bulge=True,
                 bin_array=10**np.array([8.,8.75,9.5,10.25,11.])*u.Msun):

        BinningObject.__init__(self,magnetizer_run, z=z, bins=bins,
                               bin_array=bin_array)

        self.masks = [None]*len(self.bins)

        mass = None

        if stellar_mass:
            mass = magnetizer_run.get('Mstars_disk', z)
            str_type = 'stellar'

            if include_bulge:
                mass += magnetizer_run.get('Mstars_bulge', z)
                str_disk = ' '
            else:
                str_disk = ' disc '

        if gas_mass:
            gas_mass = magnetizer_run.get('Mgas_disk', z)
            if mass is None:
                mass = gas_mass
            else:
                mass += gas_mass
            str_type = 'gas'

        if stellar_mass and gas_mass:
            str_type = 'total'

        self._quantitytype = str_type + str_disk + 'mass'
        self._update_masks(mass)



