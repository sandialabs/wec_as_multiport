
import os
import pytest

import numpy as np

import wec_as_multiport as wam
import wecopttool as wot

bem_data_fname = os.path.join(os.path.dirname(__file__),
                              'wec_as_multiport.nc')

@pytest.fixture(scope="module")
def wec():
    bem_data = wot.read_netcdf(bem_data_fname)
    wec = wam.WEC(omega=bem_data['omega'].values,
                 N=12.4666,
                 Kt=6.1745,
                 Rw=0.5,
                 Lw=0,
                 Jd=2,
                 Bd=1,
                 Kd=0,
                 Zi=np.squeeze(wot.hydrodynamic_impedance(bem_data)).values,
                 Hexc=np.squeeze(bem_data['excitation_force'].values))
    return wec


def test_low_freq_wave_hydrostatics(wec):
    """Excitation force at low frequency should approach hydrostatic 
    stiffness"""

    wave = wot.waves.regular_wave(f1=wec.f1, nfreq=wec.nfreq,
                                  freq=wec.f1, amplitude=1)

    Fexc = np.abs(wec.Fexc(wave.squeeze().values))[0]

    Kh = 0.88**2 * np.pi * 1e3 * 9.81
    Fhs = 1*Kh

    assert Fexc == pytest.approx(Fhs, rel=1e-1)


def test_Zl_opt(wec):
    """Optimal load based on bi-conjugate impedance matching condition
    should agree with optimal load based on Thevenin equivalent circuit
    """

    Zl_opt_biconj = wec.Zl_opt
    Zl_opt_Thevenin = np.conj(wec.Z_Thevenin)

    assert Zl_opt_biconj == pytest.approx(Zl_opt_Thevenin)


@pytest.mark.parametrize("wave_freq", [0.2, 0.3, 0.4, 0.5, 0.6, 0.65])
def test_max_active_power(wec, wave_freq):
    """Active power when using a load defined by the bi-conjugate matching 
    condition should match the maximum power per the Thevenin equivalent 
    system"""

    wave = wot.waves.regular_wave(f1=wec.f1, nfreq=wec.nfreq,
                                  freq=wave_freq, amplitude=1)

    Fexc = wec.Fexc(wave.squeeze().values)
    max_pow = wec.max_active_power(Fexc)
    pow = wec.active_power(Fexc)

    assert pow == pytest.approx(max_pow)
