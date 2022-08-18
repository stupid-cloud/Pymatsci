import math
import numpy as np
from scipy import constants as const
from pymatsci.processing import VaspProcess
from pymatgen.core.units import Unit, FloatWithUnit
from pymatsci.printout import print_properties, print_author_info
pi = const.pi
k = const.k
k = FloatWithUnit(k, Unit('J K^-1'))  # 玻尔兹曼常量（eV/K）
h = const.h
h = FloatWithUnit(h, Unit('J s'))  # 普朗克常量 (J/s)


class ThermalCorrection(object):
    """thermodynamic correction"""

    def __init__(self, T, p, is_line, spin_muti):
        self.T = FloatWithUnit(T, Unit('K'))  # 温度(K)
        self.p = FloatWithUnit(p, Unit('Pa'))  # 压力(Pa)
        vasp = VaspProcess()
        self.E_DFT = vasp.E_DFT
        self.sch_symbol = vasp.sch_symbol
        self.atomic_mass = vasp.atomic_mass
        self.total_atomic_mass = vasp.total_atomic_mass
        self.moment_of_inertia = vasp.moment_of_inertia
        self.sigma = vasp.sigma
        self.vib_freq = vasp.vib_freq
        self.w = spin_muti
        self.coords = vasp.coords
        self.V = k * self.T / self.p   # 气体体积
        self.is_line = is_line

    def qt(self):
        """translational partition function"""
        m = self.total_atomic_mass
        return float(math.pow(2*pi*m*k*self.T, 3/2)/math.pow(h, 3)*self.V)

    def qr1(self):
        """rotational partition function (linear molecule)"""
        I = self.moment_of_inertia_deal()
        theta = math.pow(h, 2)/(8*math.pow(pi, 2)*I*k)
        qr = 1/self.sigma*self.T/theta
        return float(qr)

    def qr2(self):
        """rotational partition function (non-linear molecule)"""
        I = self.moment_of_inertia_deal()
        theta = 1
        for i in I:
             theta *= math.pow(h, 2) / (8 * math.pow(pi, 2) * i * k)
        qr = math.pow(pi, 1/2)*math.pow(self.T, 3/2)/self.sigma/math.pow(theta, 1/2)
        return float(qr)

    def qe(self):
        """electron partition function"""
        return self.w

    def Ucorrection(self):
        """internal energy correction"""
        Et = 3/2*k*self.T    # 平动贡献
        if self.is_line:
            Er = k*self.T     # 转动贡献(线性分子)
        else:
            Er = 3/2*k*self.T  # 转动贡献(非线性分子)
        freq = self.freq_deal()

        Ev = 0    # 振动贡献
        for v in freq:  # 频率求和
            theta = h*v/k
            Ev += k*theta*(1/2+1/(math.pow(math.e, theta/self.T)-1))
        # Ev = FloatWithUnit(Ev, Unit('J'))
        return FloatWithUnit(Et+Er+Ev, Unit('J'))

    def Scorrection(self):
        """entropy correction"""
        St = k*(math.log(self.qt())+3/2)    # 平动贡献
        if self.is_line:
            Sr = k*(math.log(self.qr1())+1)     # 转动贡献(线性分子)
        else:
            Sr = k * (math.log(self.qr2())+3/2)  # 转动贡献(非线性分子)
        Sv = 0      # 振动贡献
        freq = self.freq_deal()
        for v in freq:  # 频率求和
            theta = h*v/k
            Sv += theta/self.T/(math.pow(math.e, theta/self.T)-1)-math.log(1-math.pow(math.e, -theta/self.T))
        Se = k*(math.log(self.qe())+0)     # 电子贡献
        return St+Sr+Se+k

    def Zcorrection(self):
        """zero point energy"""
        m = self.total_atomic_mass

        m = float(m)
        Zt =3*h*h/8/m/math.pow(self.V, 2/3)  # 平动零点能
        freq = self.freq_deal()
        fre = 0  # 振动零点能
        for v in freq:  # 频率求和
            fre += v
        Zv = 1/2*fre*h
        # print(Zt .to('eV'), Zv.to('eV'))
        return FloatWithUnit(Zv+Zt, Unit('J'))

    def freq_deal(self):
        """Imaginary frequency and related processing of vibration degrees of freedom"""
        real_freq, imag_freq, real_num, imag_num = self.vib_freq  # 计算的得到的频率
        real_freq.sort()
        freq = real_freq  # 虚频太多时
        if self.is_line:
            if imag_num < 5:
                freq = real_freq[5 - imag_num:]
        else:
            if imag_num < 6:
                freq = real_freq[6 - imag_num:]
        return freq

    def moment_of_inertia_deal(self):
        """processing of moment of inertia"""

        eig_val, eig_vec = np.linalg.eigh(self.moment_of_inertia)
        d = np.dot(np.dot(np.linalg.inv(eig_vec), self.moment_of_inertia), eig_vec)
        # 根据是否是线性分子给出I
        I_x = d[0, 0]
        I_y = d[1, 1]
        I_z = d[2, 2]

        if self.is_line:
            I = (I_x+I_y+I_z)/2
        else:
            I = [I_x, I_y, I_z]
        return I

    def printout(self):
        """print information"""
        properties1 = {"Temperature (T)": self.T, "Pressure (P)": self.p, "Linear molecule": self.is_line, "Spin-Multiplicity (S)": self.w,
                      "Point group": self.sch_symbol, 'symmetry number': self.sigma}
        properties2 = {'Energy of DFT (E_DFT)': self.E_DFT, "Zero-point energy (E_ZPE)": self.Uz,
                      "Thermal correction to U(T)": self.U, "Thermal correction to H(T)": self.H, "Thermal correction to G(T)": self.G, "Entropy S": self.S,
                      "Entropy contribution T*S": self.T * self.S}
        properties3 = {"corrected E_DFT": self.E_DFT+self.Uz, "Corrected U(T)": self.E_DFT+self.Uz+self.U,
                      "Corrected H(T)": self.E_DFT+self.Uz+self.H, "Corrected G(T)": self.E_DFT+self.Uz+self.G}
        print_author_info()
        print_properties(properties1)
        print_properties(properties2)
        print_properties(properties3)


    def free_gas_correction(self):
        """free gas correction"""
        self.Uz = self.Zcorrection().to('eV')
        self.S = self.Scorrection().to('eV K^-1')
        self.U = self.Ucorrection().to('eV')
        self.H = self.U + (self.p*self.V).to('eV')
        self.G = self.H - self.T * self.S



if __name__ == '__main__':

    # 选择自由气体或者吸附气体的修正
    # print("The thermodynamic correction includes two cases: 1. Free gas; 2. Gas adsorbed on the surface")
    # slelect = input('Please enter your choice')

    corre = ThermalCorrection()
    # corre.Zcorrection()
    corre.correction()  # 进行修正
    corre.printout()    # 打印输出

