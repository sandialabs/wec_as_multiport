# © 2024 National Technology & Engineering Solutions of Sandia, LLC
# (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
# Government retains certain rights in this software.

import numpy as np
import copy

from wec_as_multiport import util

# import capytaine as cpy
# import wecopttool as wot

__all__ = [
    "WEC",
    "__max_active_power__",
    "__F_Thevenin__",
    "__power__",
    "__active_power__",
    "__reactive_power__",
    "__apparent_power__",
    "__Fexc__",
    "__Zin__",
    "__Zout__",
    "__pid_controller__",
    "power_reflection_coefficient",
    "power_transmission_coefficient",
]


class WEC:

    def __init__(self, omega, N, Kt, Rw, Lw, Jd, Bd, Kd, Zi, Hexc, name=None) -> None:

        self.omega = omega
        self.N = N
        self.Kt = Kt
        self.Rw = Rw
        self.Lw = Lw
        self.Jd = Jd
        self.Bd = Bd
        self.Kd = Kd

        self.Zi = Zi
        self.Hexc = Hexc
        
        if name is None:
            name = ''
        self.name = name

    def __repr__(self) -> str:
        type_ = type(self)
        module = type_.__module__
        qualname = type_.__qualname__
        repr_org = f"<{module}.{qualname} object at {hex(id(self))}>"
        rs = repr_org + "\n" + \
            f'\t\t{self.name}\n' + \
            f'\tHydro. resonance:\t{self.hydrodynamic_resonance:.2f}Hz\n' + \
            f'\tThèvenin resonance:\t{self.Thevenin_resonance:.2f}Hz\n'
        return rs

    @property
    def f1(self) -> float:
        return self.omega[0]/2/np.pi

    @property
    def nfreq(self) -> int:
        """Number of frequencies"""
        return len(self.omega)

    @property
    def freq(self) -> np.ndarray:
        return self.omega/2/np.pi

    @property
    def Zd(self) -> np.ndarray:
        """Drive train impedance"""
        return 1j * self.omega * self.Jd + self.Bd \
            + self.Kd / (1j*self.omega)

    @property
    def Zw(self) -> np.ndarray:
        """Winding impedance"""
        return self.Rw + 1j * self.omega * self.Lw

    @property
    def Zpto(self) -> np.ndarray:
        """PTO matrix in impedance form"""
        return np.array([
            [self.Zd*self.N**2, -self.Kt*self.N*np.ones_like(self.omega)],
            [self.Kt*self.N*np.ones_like(self.omega), self.Zw]])

    @property
    def Ypto(self) -> np.ndarray:
        """PTO matrix in admittance form"""
        return 1/self.__detZpto__ * np.array([
            [self.Zpto[1, 1], -1*self.Zpto[0, 1]],
            [-1*self.Zpto[1, 0], self.Zpto[0, 0]]])

    @property
    def ABCDpto(self) -> np.ndarray:
        """PTO matrix in ABCD form"""
        return 1/self.Zpto[1, 0] * np.array([
            [self.Zpto[0, 0], self.__detZpto__],
            [np.ones_like(self.omega), self.Zpto[1, 1]]])

    @property
    def Hpto(self) -> np.ndarray:
        """PTO matrix in hybrid form"""
        return 1/self.Zpto[1, 1] * np.array([
            [self.__detZpto__, self.Zpto[0, 1]],
            [-1*self.Zpto[1, 0], np.ones_like(self.omega)]])

    @property
    def Zout(self) -> np.ndarray:
        """PTO output impedance"""
        return __Zout__(self.Zpto, self.Zi)

    @property
    def Zl_opt(self) -> np.ndarray:
        """Optimal load impedance"""
        return np.conj(self.Zout)

    @property
    def hydrodynamic_resonance_index(self) -> int:
        """Index of hydrodynamic resonance frequency"""
        return np.argmin(np.abs(np.angle(self.Zi)))

    @property
    def hydrodynamic_resonance(self) -> float:
        """Hydrodynamic resonant frequency"""
        # return self.freq[self.hydrodynamic_resonance_index]
        try:
            fn = util.find_zero_crossings(self.freq, np.angle(self.Zi))[0]
        except:
            fn = self.freq[self.hydrodynamic_resonance_index]
            # TODO add warning

        return fn

    @property
    def Zl_opt_mech(self) -> np.ndarray:
        """Load impedance for optimal mechanical power"""
        Zlm = np.conj(self.Zi)
        return self.Zpto[0, 1]*self.Zpto[1, 0] / (self.Zpto[0, 0] - Zlm) \
            - self.Zpto[1, 1]

    @property
    def Z_Thevenin(self) -> np.ndarray:
        """Thévenin impedance"""
        return self.Zout

    @property
    def Thevenin_resonance_index(self) -> int:
        """Index of Thévenin resonance frequency"""
        return np.argmin(np.abs(np.angle(self.Z_Thevenin)))

    @property
    def Thevenin_resonance(self) -> float:
        """Thévenin resonant frequency"""
        # return self.freq[self.Thevenin_resonance_index]
        return util.find_zero_crossings(self.freq, np.angle(self.Z_Thevenin))[0]

    @property
    def __detZpto__(self) -> np.array:
        """Determinant of Zpto"""
        return np.linalg.det(np.transpose(self.Zpto, [2, 0, 1]))

    @property
    def Hexc_Thevenin(self) -> np.array:
        """Thévenin excitation transfer function"""
        return __H_Thevenin__(self.Zpto, self.Zi, self.Hexc)

    # def Spto(self, Zl=None): -> np.ndarray:
    #     "Power scattering matrix per TODO"

    #     if Zl is None:
    #         Zl = self.Zl_opt

    #     S11 = (self.Zpto[0,0] - np.conj(self.Zi))*(self.Zpto[1,1] + Zl) \
    #         - (self.Zpto[0,1]*self.Zpto[1,0])

    def transducer_power_gain(self, Zl=None) -> np.ndarray:
        """Wave-to-wire efficiency (input power / power delivered to load)
        https://en.wikipedia.org/wiki/Power_gain#Transducer_power_gain"""
        if Zl is None:
            Zl = self.Zl_opt
        return 4*np.abs(self.Zpto[1, 0])**2 * np.real(Zl)*np.real(self.Zi) \
            / np.abs((self.Zpto[0, 0] + self.Zi)*(self.Zpto[1, 1] + Zl)
                     - self.Zpto[0, 1]*self.Zpto[1, 0])**2

    def available_power_gain(self) -> np.ndarray:
        """Ratio of maximum power delivered to the load to the maximum power 
        available at the source (i.e., the “ideal wave-to-wire efficiency”)"""
        return np.abs(self.Zpto[1, 0]/(self.Zi + self.Zpto[0, 0]))**2 * \
            np.real(self.Zi) / np.real(self.Zout)

    def operating_power_gain(self, Zl=None) -> np.ndarray:
        """Ratio of power delivered to the load to the power delivered to the 
        power take-off (i.e., the “PTO efficiency”)"""
        if Zl is None:
            Zl = self.Zl_opt
        return np.abs(self.Zpto[1, 0]/(Zl + self.Zpto[1, 1]))**2 * \
            np.real(Zl) / np.real(self.Zin(Zl=Zl))

    def F_Thevenin(self, Fexc) -> np.ndarray:
        """Thévenin source"""
        return __F_Thevenin__(self.Zpto, self.Zi, Fexc)

    def Zin(self, Zl=None) -> np.ndarray:
        """PTO input impedance"""
        if Zl is None:
            Zl = self.Zl_opt
        return __Zin__(self.Zpto, Zl)

    def Fexc(self, waves) -> np.ndarray:
        """Excitation spectrum"""
        return __Fexc__(self.Hexc, waves)

    def power_variables_in(self, Fexc, Zl=None) -> np.ndarray:
        """Power variables before Zpto"""
        if Zl is None:
            Zl = self.Zl_opt
        # Fpto = self.Fpto(Fexc=Fexc, Zl=Zl)
        # v = self.velocity(Fexc=Fexc, Zl=Zl)
        v = Fexc / (self.Zi + self.Zin(Zl=Zl))
        Fpto = v * self.Zin(Zl=Zl)
        return np.array([[Fpto], [v]])

    def power_variables_out(self, Fexc, Zl=None):
        """Power variables after Zpto"""
        if Zl is None:
            Zl = self.Zl_opt
        vars_in = self.power_variables_in(Fexc, Zl=Zl)
        # return np.einsum('mnf,nkf->mkf', self.invABCDpto, vars_in)
        return np.einsum('mnf,nkf->mkf', self.Bpto, vars_in)

    def power_mech(self, Fexc, Zl=None):
        """Complex power at PTO input"""
        Fpto, v = self.power_variables_in(Fexc=Fexc, Zl=Zl)
        return 1/2*np.conj(np.squeeze(Fpto))*np.squeeze(v)

    def active_power_mech(self, Fexc, Zl=None):
        """Active power at PTO input"""
        return __active_power__(self.power_mech(Fexc, Zl))

    def reactive_power_mech(self, Fexc, Zl=None):
        """Reactive power at PTO input"""
        return __reactive_power__(self.power_mech(Fexc, Zl))

    def apparent_power_mech(self, Fexc, Zl=None):
        """Apparent power at PTO input"""
        return __apparent_power__(self.power_mech(Fexc, Zl))

    def power(self, Fexc, Zl=None):
        """Complex power at load"""
        if Zl is None:
            Zl = self.Zl_opt
        return __power__(self.Zpto, self.Zi, Fexc, Zl)

    def active_power(self, Fexc, Zl=None):
        """Active power at load"""
        return __active_power__(self.power(Fexc, Zl))

    def reactive_power(self, Fexc, Zl=None):
        """Reactive power at load"""
        return __reactive_power__(self.power(Fexc, Zl))

    def apparent_power(self, Fexc, Zl=None):
        """Apparent power at load"""
        return __apparent_power__(self.power(Fexc, Zl))

    def max_active_power(self, Fexc):
        """Maximum active power"""
        return __max_active_power__(self.Z_Thevenin, self.F_Thevenin(Fexc))

    def max_active_power_mech(self, Fexc):
        """Maximum active mechanical power"""
        return __max_active_power__(self.Zi, Fexc)

    def pi_opt(self, freq) -> tuple:
        """Optimal PI gains for a controller acting on current and shaft speed
        for a given frequency"""
        indx = np.argmin(np.abs(self.freq - freq))
        Zc_opt = self.Kt/(np.conj(self.Zout) + self.Zw)
        kp = np.real(Zc_opt[indx])
        ki = np.real(Zc_opt[indx]*1j*self.omega[indx])
        return (kp, ki)

    def pid_controller(self, kp=0, ki=0, kd=0) -> np.ndarray:
        """Controller impedance"""
        return __pid_controller__(self.omega, kp, ki, kd)

    def Zl_C(self, C) -> np.ndarray:
        """Load impedance due to feedback controller acting on current and 
        shaft speed"""
        return self.Kt/C - self.Zw

    def copy(self):
        return copy.deepcopy(self)


def __max_active_power__(Z, F):
    return np.abs(F)**2 / (8 * np.real(Z))


def __H_Thevenin__(Zpto, Zi, Hexc):
    return Hexc*Zpto[1, 0] / (Zi + Zpto[0, 0])


def __F_Thevenin__(Zpto, Zi, Fexc):
    return Fexc*Zpto[1, 0] / (Zi + Zpto[0, 0])


def __power__(Zpto, Zi, Fexc, Zl=None) -> np.ndarray:
    # TODO - define in terms of flow and effort
    return np.abs(Fexc*Zpto[1, 0]
                  / ((Zpto[1, 1] + Zl)*(Zi + Zpto[0, 0]) -
                      Zpto[1, 0]*Zpto[0, 1]))**2 * 1/2*Zl


def __active_power__(power):
    return np.real(power)


def __reactive_power__(power):
    return np.im(power)


def __apparent_power__(power):
    return np.abs(power)


def __Fexc__(Hexc, waves) -> np.ndarray:
    return Hexc * waves


def __Zin__(Z_2port, Zl) -> np.ndarray:
    return Z_2port[0, 0] - Z_2port[0, 1] * Z_2port[1, 0] \
        / (Z_2port[1, 1] + Zl)


def __Zout__(Z_2port, Zi) -> np.ndarray:
    return Z_2port[1, 1] - Z_2port[1, 0] * Z_2port[0, 1] \
        / (Zi + Z_2port[0, 0])


def __pid_controller__(omega, kp=0, ki=0, kd=0) -> np.ndarray:
    return kp + ki/(1j*omega) + kd*1j*omega


def power_reflection_coefficient(Zs, Zl) -> np.ndarray:
    """Power reflection coefficient per Kurokawa 1965 eq. 14

    K. Kurokawa, "Power Waves and the Scattering Matrix," in IEEE Transactions 
    on Microwave Theory and Techniques, vol. 13, no. 2, pp. 194-202, March 
    1965, doi: 10.1109/TMTT.1965.1125964.
    """
    s = (Zl - np.conj(Zs))/(Zl + Zs)
    return np.abs(s)**2


def power_transmission_coefficient(Zs, Zl) -> np.ndarray:
    """Power transmission coefficient per Kurokawa 1965

    K. Kurokawa, "Power Waves and the Scattering Matrix," in IEEE Transactions 
    on Microwave Theory and Techniques, vol. 13, no. 2, pp. 194-202, March 
    1965, doi: 10.1109/TMTT.1965.1125964.
    """

    return 1 - power_reflection_coefficient(Zs=Zs, Zl=Zl)
