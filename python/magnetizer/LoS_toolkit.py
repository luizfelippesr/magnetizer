import numpy as np
import subprocess
import os

def run_LoS(runtype, parameters_file, gal_id, iz, theta, wavelength=20e-2,
            ignore_small_scale=False, zmax=None, ymax=None, imgdir='img',
            cwd=None, verbose=False):
    theta_rad = np.deg2rad(theta)

    os.makedirs(imgdir, exist_ok=True)

    cmd = "./Observables_single.exe {type} {param} {ig} {iz} {theta} {lamb} {ignore_small}".format(
      param=parameters_file, ig=gal_id+1, iz=iz, theta=theta_rad,
      lamb=wavelength, ignore_small=int(ignore_small_scale), type=runtype)

    if zmax is not None:
        zmax_actual = zmax/np.sin(theta_rad)
        cmd += " {ymax} {zmax} {dir}".format(ymax=ymax, zmax=zmax_actual,
                                             dir=imgdir)

    if 'RM_study' in runtype:
        cmd += ' > /tmp/RMstudy.tmp'

    cmd = cmd.split(' ')
    if verbose:
        print(' '.join(cmd))
    try:
        output = subprocess.check_output(cmd, cwd=cwd)
    except subprocess.CalledProcessError as e:
        print('error')
        print(' '.join(cmd))
        print(e.output.decode())
        print('---')
        raise

    if verbose:
        print( output.decode() )
    if 'RM_study' in runtype:
        output = np.genfromtxt('/tmp/RMstudy.tmp',skip_header=3, unpack=True)

    return output


class Stokes_data(object):
    def __init__(self, parameters_file, gal_id, iz, theta, wavelength=20e-2,
                 ignore_small_scale=False, zmax=None, dust=False, test=False,
                 ymax=None, imgdir='/tmp/img', cwd=None, verbose=False):
        if test:
            imageType = 'testImage'
        elif not dust:
            imageType = 'Image'
        else:
            imageType = 'dustImage'

        run_LoS(imageType, parameters_file, gal_id, iz, theta, wavelength,
                ignore_small_scale=ignore_small_scale,
                zmax=zmax, ymax=ymax, cwd=cwd,
                imgdir=imgdir, verbose=verbose)

        for d in ['z','y','I','Q','U','RM','cells']:
            path = os.path.join(imgdir, d+'.dat')
            data = np.genfromtxt(path)
            setattr(self, d, np.genfromtxt(path))
        self.theta = theta
        self.wavelength = wavelength
        self._PI = None
        self._Psi = None
        self._p = None
        self._zsky = None
        
        self.extent = (self.y.min(), self.y.max(), 
                       self.z_sky.min(), self.z_sky.max())
        
    @property
    def z_sky(self):
        if self._zsky is None:
            self._zsky = self.z*np.sin(np.deg2rad(self.theta))
        return self._zsky

    @property
    def PI(self):
        if self._PI is None:
            self._PI = np.sqrt(self.Q**2+self.U**2)
        return self._PI
    @property
    def p(self):
        if self._p is None:
            self._p = np.abs(self.PI/self.I)
        return self._p

    @property
    def Psi(self):
        if self._Psi is None:
            self._Psi = 0.5*np.arctan2(self.U,self.Q)
        return self._Psi


class RM_study(object):
    def __init__(self, parameters_file, gal_id, iz, theta,
                 ignore_small_scale=False, verbose=False):
        from StringIO import StringIO
        self.RM, self.r, self.y, self.z = run_LoS('RM_study', parameters_file,
                                                  gal_id, iz, theta,
                                                  ignore_small_scale=ignore_small_scale,
                                                  verbose=verbose)
