
import numpy as np
import copy

# import capytaine as cpy
# import wecopttool as wot


class WEC:

    def __init__(self, omega, N, Kt, Rw, Lw, Jd, Bd, Kd, Zi, Hexc) -> None:

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

    # def __repr__(self) -> str:
    #     type_ = type(self)
    #     module = type_.__module__
    #     qualname = type_.__qualname__
    #     repr_org = f"<{module}.{qualname} object at {hex(id(self))}>"
    #     return repr_org + " :: " + self.__str__()

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
        return np.argmin(np.abs(np.imag(self.Zi)))

    @property
    def hydrodynamic_resonance(self) -> float:
        """Hydrodynamic resonant frequency"""
        return self.freq[self.hydrodynamic_resonance_index]

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
    def __detZpto__(self) -> np.array:
        """Determinant of Zpto"""
        return np.linalg.det(np.transpose(self.Zpto, [2, 0, 1]))

    # TODO
    # def Fpto(self, Fexc, Zl=None) -> np.ndarray:
    #     """PTO force"""
    #     if Zl is None:
    #         Fpto = self.velocity(Fexc=Fexc, Fpto=None) * self.Zin(Zl=Zl)
    #     else:
    #         raise NotImplementedError() #TODO
    #     return Fpto

    # def velocity(self, Fexc, Fpto=None) -> np.ndarray:
    #     """Rectilinear velocity"""
    #     if Fpto is None:
    #         v = Fexc / (self.Zi + self.Zin(Zl=self.Zl_opt))
    #     else:
    #         v = (Fexc - Fpto) / self.Zi
    #     return v

    def transducer_power_gain(self, Zl=None) -> np.ndarray:
        """Wave-to-wire efficiency (input power / power delivered to load)
        https://en.wikipedia.org/wiki/Power_gain#Transducer_power_gain"""
        if Zl is None:
            Zl = self.Zl_opt
        return 4*np.abs(self.Zpto[1, 0])**2 * np.real(Zl)*np.real(self.Zi) \
            / np.abs((self.Zpto[0, 0] + self.Zi)*(self.Zpto[1, 1] + Zl) \
                - self.Zpto[0, 1]*self.Zpto[1, 0])**2

    def Zlm(self, Zl) -> np.ndarray:
        """Mechanical load impedance"""
        return self.Zin(Zl=Zl)

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

    def power_mech(self, Fexc, Zl=None):
        if Zl is None:
            Zl = self.Zl_opt
        raise NotImplementedError()  # TODO

    def power(self, Fexc, Zl=None):
        """Complex power at load"""
        # TODO - define in terms of Vout and Iout
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

    def copy(self):
        return copy.deepcopy(self)


def __max_active_power__(Z, F):
    return np.abs(F)**2 / (8 * np.real(Z))


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


def __Zin__(Zpto, Zl) -> np.ndarray:
    return Zpto[0, 0] - Zpto[0, 1] * Zpto[1, 0] / (Zpto[1, 1] + Zl)


def __Zout__(Zpto, Zi) -> np.ndarray:
    return Zpto[1, 1] - Zpto[1, 0] * Zpto[0, 1] / (Zi + Zpto[0, 0])
