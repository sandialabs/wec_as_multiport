
import os
import pytest

import numpy as np

import wec_as_multiport as wam
import wecopttool as wot
import capytaine as cpy

bem_data_fname = os.path.join(os.path.dirname(__file__),
                              '..', 'data', 'wec_as_multiport.nc')


@pytest.fixture(scope="module")
def wec(Rw=None):
    if Rw is None:
        Rw = 0.5

    if os.path.isfile(bem_data_fname):
        print("Found existing BEM file, loading")
        bem_data = wot.read_netcdf(bem_data_fname)
    else:
        f1 = 0.025
        nfreq = 60
        freq = wot.frequency(f1, nfreq, False)  # False -> no zero frequency

        wb = wot.geom.WaveBot()  # use standard dimensions
        mesh_size_factor = 0.5  # 1.0 for default, smaller to refine mesh
        mesh = wb.mesh(mesh_size_factor)
        fb = cpy.FloatingBody.from_meshio(mesh, name="WaveBot")
        fb.add_translation_dof(name="Heave")
        bem_data = wot.run_bem(fb, freq)
        bem_data = bem_data.assign_coords(
            freq=("omega", bem_data['omega'].values/2/np.pi))
        bem_data['freq'].attrs['long_name'] = 'Frequency'
        bem_data['freq'].attrs['units'] = 'Hz'
        bem_data['excitation_force'] = bem_data['diffraction_force'] + \
            bem_data['Froude_Krylov_force']
        bem_data = wot.add_linear_friction(bem_data)

    wec = wam.WEC(omega=bem_data['omega'].values,
                  N=12.4666,
                  Kt=6.1745,
                  Rw=Rw,
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


def test_Zlm_matching(wec):
    """Mechanical load impedance when maximizing mechanical power should 
    equal the complex conjugate of the intrinsic impedance"""

    Z1 = wec.Zlm(Zl=wec.Zl_opt_mech)
    Z2 = np.conj(wec.Zi)

    np.testing.assert_allclose(Z1, Z2)


def test_Zl_opt(wec):
    """Optimal load based on bi-conjugate impedance matching condition
    should agree with optimal load based on Thevenin equivalent circuit
    """

    Zl_opt_biconj = wec.Zl_opt
    Zl_opt_Thevenin = np.conj(wec.Z_Thevenin)

    np.testing.assert_allclose(Zl_opt_biconj, Zl_opt_Thevenin)


def test_no_reflection(wec):
    """Expect no power reflection with perfect matching"""

    Gamma = wam.power_reflection_coefficient(Zs=wec.Z_Thevenin,
                                             Zl=wec.Zl_opt)

    np.testing.assert_allclose(np.zeros_like(Gamma), Gamma)


def test_reflection(wec):
    """Should see power reflection with imperfect matching"""
    
    Gamma = wam.power_reflection_coefficient(Zs=wec.Z_Thevenin,
                                             Zl=wec.Zl_opt_mech)

    np.testing.assert_array_less(np.zeros_like(Gamma), Gamma)


@pytest.mark.parametrize("design_freq", [0.4, 0.55, 0.65])
def test_reflection_in_bounds(wec,design_freq):
    """This power reflection coefficient must be between 0 and 1"""

    kp, ki = wec.pi_opt(design_freq)
    C = wec.pid_controller(kp=kp, ki=ki)
    Zl_C = wec.Zl_C(C)
    gamma_pi = wam.power_reflection_coefficient(Zs=wec.Zi,
                                                Zl=wec.Zin(Zl=Zl_C))

    np.testing.assert_array_less(np.zeros_like(gamma_pi), gamma_pi)
    np.testing.assert_array_less(gamma_pi, np.ones_like(gamma_pi))

@pytest.mark.parametrize("wec", np.logspace(-5, -1, 5), indirect=True)
class TestPowerGains:

    def test_transducer_equal_available_power_gain(self, wec):
        """Transducer gain should match available power gain if optimal load is 
        used"""

        apg = wec.available_power_gain()
        tpg = wec.transducer_power_gain(Zl=wec.Zl_opt)

        np.testing.assert_allclose(apg, tpg)

    def test_available_power_gain_less_than_unity(self, wec):
        """Available power gain should always be less than 1 if you have losses"""

        apg = wec.available_power_gain()

        assert np.all(apg < 1)

    def test_transducer_equal_available_power_gain(self, wec):
        """Transducer gain should match available power gain if optimal load is 
        used"""

        apg = wec.available_power_gain()
        tpg = wec.transducer_power_gain(Zl=wec.Zl_opt)

        np.testing.assert_allclose(apg, tpg)

    def test_operating_power_gain_less_than_unity(self, wec):
        """Operating power gain should always be less than unity if you have 
        losses"""

        opg = wec.operating_power_gain()

        assert np.all(opg < 1)


def test_operating_power_gain_is_unity(wec):
    """Operating power gain should be unity if there are no losses"""
    wec.Rw = 0
    wec.Bd = 0

    opg = wec.operating_power_gain()

    np.testing.assert_allclose(np.ones_like(opg), opg, rtol=1e-14)


@pytest.mark.parametrize("wave_freq", np.linspace(0.2, 0.65, 10))
class TestPerfomanceAtFreqs:
    def test_max_active_power(self, wec, wave_freq):
        """Active power when using a load defined by the bi-conjugate matching 
        condition should match the maximum power per the Thevenin equivalent 
        system"""

        wave = wot.waves.regular_wave(f1=wec.f1, nfreq=wec.nfreq,
                                      freq=wave_freq, amplitude=1)

        Fexc = wec.Fexc(wave.squeeze().values)
        max_pow = wec.max_active_power(Fexc)
        pow = wec.active_power(Fexc)

        assert pow == pytest.approx(max_pow)

    @pytest.fixture()
    def freq_ind(self, wec, wave_freq):
        return np.argmin(np.abs(wec.freq - wave_freq))

    @pytest.fixture()
    def pi_load_impedance(self, wec, wave_freq):
        kp, ki = wec.pi_opt(freq=wave_freq)
        C = wec.pid_controller(kp=kp, ki=ki)
        return wec.Zl_C(C)

    def test_pi_load_matches_optimal_load_at_design_freq(self, wec,
                                                         freq_ind,
                                                         pi_load_impedance):

        assert pi_load_impedance[freq_ind] == pytest.approx(
            wec.Zl_opt[freq_ind])

    def test_pi_transducer_power_gain_matches_optimal_at_design_freq(sef,
                                                                     wec,
                                                                     freq_ind,
                                                                     pi_load_impedance):

        tpg = wec.transducer_power_gain(Zl=pi_load_impedance)
        apg = wec.available_power_gain()
        assert tpg[freq_ind] == pytest.approx(apg[freq_ind])
