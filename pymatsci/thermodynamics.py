import math
import numpy as np
from scipy import constants as const
from pymatsci.processing import VaspProcess
from pymatgen.core.units import Unit, FloatWithUnit
from pymatsci.printout import print_properties, print_author_info
import requests
import re


pi = const.pi
k = const.k
k = FloatWithUnit(k, Unit('J K^-1'))  # 玻尔兹曼常量（eV/K）
h = const.h
h = FloatWithUnit(h, Unit('J s'))  # 普朗克常量 (J/s)
Na = const.Avogadro


class FreeGasCorrection(object):
    """free gas correction"""

    def __init__(self, T, p, is_line, spin_muti):
        self.T = FloatWithUnit(T, Unit('K'))  # 温度(K)
        self.p = FloatWithUnit(p, Unit('Pa'))  # 压力(Pa)
        vasp = VaspProcess()
        self.E_DFT = vasp.E_DFT
        self.sch_symbol = vasp.sch_symbol
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
        Ev = 0
        for v in freq:  # 频率求和
            theta = h*v/k
            Ev += k*theta*(1/(math.pow(math.e, theta/self.T)-1))  # 振动贡献
        # Ev = FloatWithUnit(Ev, Unit('J'))
        Ee = 0  # 电子贡献
        return FloatWithUnit(Et+Er+Ev+Ee, Unit('J'))

    def Scorrection(self):
        """entropy correction"""
        St = k*(math.log(self.qt())+3/2)    # 平动贡献
        if self.is_line:
            Sr = k*(math.log(self.qr1())+1)     # 转动贡献(线性分子)
        else:
            Sr = k * (math.log(self.qr2())+3/2)  # 转动贡献(非线性分子)

        Sv = 0
        freq = self.freq_deal()
        for v in freq:  # 频率求和
            theta = h*v/k
            Sv += theta/self.T/(math.pow(math.e, theta/self.T)-1)-math.log(1-math.pow(math.e, -theta/self.T))
        Sv = Sv * k   # 振动贡献
        Se = k*(math.log(self.qe())+0)     # 电子贡献
        return St+Sr+Se+Sv+k

    def Zcorrection(self):
        """zero point energy"""

        freq = self.freq_deal()
        fre = 0  # 振动零点能
        for v in freq:  # 频率求和
            fre += v
        Zv = 1/2*fre*h
        # print(Zt .to('eV'), Zv.to('eV'))
        return FloatWithUnit(Zv, Unit('J'))

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
        properties3 = {"corrected E_DFT": self.E_DFT+self.Uz, "Corrected U(T)": self.E_DFT+self.U,
                      "Corrected H(T)": self.E_DFT+self.H, "Corrected G(T)": self.E_DFT+self.G}
        print_author_info()
        print_properties(properties1)
        print_properties(properties2)
        print_properties(properties3)

    def correction(self):
        """free gas correction"""
        self.Uz = self.Zcorrection().to('eV')
        self.S = self.Scorrection().to('eV K^-1')
        self.U = self.Ucorrection().to('eV') + self.Uz
        self.H = self.U + (self.p*self.V).to('eV')
        self.G = self.H - self.T * self.S


class AdsorbedGasCorrection(object):
    """adsorbed gas correction"""

    def __init__(self, T):
        self.T = FloatWithUnit(T, Unit('K'))  # 温度(K)
        vasp = VaspProcess()
        self.vib_freq = vasp.vib_freq
        self.E_DFT = vasp.E_DFT

    def Ucorrection(self):
        """internal energy correction"""
        freq = self.freq_deal()
        Ev = 0    # 振动贡献
        for v in freq:  # 频率求和
            theta = h*v/k
            Ev += k*theta*(1/(math.pow(math.e, theta/self.T)-1))
        # Ev = FloatWithUnit(Ev, Unit('J'))
        return FloatWithUnit(Ev, Unit('J'))

    def Scorrection(self):
        """entropy correction"""
        # 振动贡献
        Sv = 0
        freq = self.freq_deal()
        for v in freq:  # 频率求和
            theta = h*v/k
            Sv += theta/self.T/(math.pow(math.e, theta/self.T)-1)-math.log(1-math.pow(math.e, -theta/self.T))
        Sv = k * Sv
        return FloatWithUnit(Sv, Unit('J K^-1'))

    def Zcorrection(self):
        """zero point energy"""

        freq = self.freq_deal()
        fre = 0  # 振动零点能
        for v in freq:  # 频率求和
            fre += v
        Zv = 1/2*fre*h
        # print(Zt .to('eV'), Zv.to('eV'))
        return FloatWithUnit(Zv, Unit('J'))

    def freq_deal(self):
        """Imaginary frequency and related processing of vibration degrees of freedom"""
        real_freq, imag_freq, real_num, imag_num = self.vib_freq  # 计算的得到的频率
        for index in range(len(real_freq)):
            if real_freq[index] < 1.49896E12:
                real_freq[index] = 1.49896E12

        return real_freq

    def printout(self):
        """print information"""
        properties1 = {"Temperature (T)": self.T}
        properties2 = {'Energy of DFT (E_DFT)': self.E_DFT, "Zero-point energy (E_ZPE)": self.Uz,
                      "Thermal correction to U(T)": self.U, "Thermal correction to H(T)": self.H, "Thermal correction to G(T)": self.G, "Entropy S": self.S,
                      "Entropy contribution T*S": self.T * self.S}
        properties3 = {"corrected E_DFT": self.E_DFT+self.Uz, "Corrected U(T)": self.E_DFT+self.U,
                      "Corrected H(T)": self.E_DFT+self.H, "Corrected G(T)": self.E_DFT+self.G}
        print_author_info()
        print_properties(properties1)
        print_properties(properties2)
        print_properties(properties3)

    def correction(self):
        """free gas correction"""
        self.Uz = self.Zcorrection().to('eV')
        self.S = self.Scorrection().to('eV K^-1')
        self.U = self.Ucorrection().to('eV') + self.Uz
        self.H = self.U
        self.G = self.H - self.T * self.S


class Abthermodynamics(object):
    """从头算热力学"""
    def __init__(self, url):
        self.url = url

    def therm_table(self, url):
        """获取热力学表"""
        r = requests.get(url)
        # print(r.text)
        data = np.array(re.findall('<TD>(.+)</TD>', r.text))
        data = data.reshape((int(len(data) / 8)), 8)
        return data

    def to_dict(self, url):
        data = self.therm_table(url)
        dict_data = dict()
        for i in range(len(data[:, 0])):
            dict_data[data[:, 0][i]] = data[i, 1:]
        return dict_data

    def get_T(self, n, E_ads, p):
        """固定p，得到T"""
        table = self.therm_table(url=self.url)
        data = table[:, [0, 4]].astype('float')
        T = table[:, 0].astype('float')
        S = table[:, 2].astype('float')
        HT_H = table[:, 4].astype('float')
        A = E_ads-n*(HT_H-HT_H[0]-T*S/1000)/Na*1000/1.602177E-19
        B = float(n*k.to('eV K^-1')*math.log(p))*T
        data[:, 0] = T
        data[:, 1] = A-B
        return data

    def get_p(self, n, E_ads, T):
        """固定T，得到p"""
        dict_data = self.to_dict(self.url)
        S = float(dict_data[str(T)][1])
        dH = float(dict_data[str(T)][3])-float(dict_data[str(0)][3])
        A = E_ads - n * (dH - T * S / 1000) / Na * 1000 / 1.602177E-19
        B = float(n * k.to('eV K^-1') * T)
        return '{}-({}*T)'.format(A, B)

    def get_Tp(self, n1, E_ads1, n2, E_ads2):
        table = self.therm_table(self.url)
        data = table[1:, [0, 4]].astype('float')
        T = table[1:, 0].astype('float')
        S = table[1:, 2].astype('float')
        HT_H = table[1:, 4].astype('float')
        A1 = E_ads1 - n1 * (HT_H - HT_H[0] - T * S / 1000) / Na * 1000 / 1.602177E-19
        A2 = E_ads2 - n2 * (HT_H - HT_H[0] - T * S / 1000) / Na * 1000 / 1.602177E-19
        B1 = -n1 * float(k.to('eV K^-1')) * T
        B2 = -n2 * float(k.to('eV K^-1')) * T
        p = (A1-A2)/(B2-B1)
        data[:, 0] = T
        data[:, 1] = p
        return data


if __name__ == '__main__':

    t = Abthermodynamics('https://janaf.nist.gov/tables/O-029.html')
    data1 = t.get_p(1/2, -5.134, 300)  # 输入覆盖度，吸附能，温度
    data2 = t.get_p(2/2, -10.652, 300)  # 输入覆盖度，吸附能，温度
    data = np.hstack((data1, data2))
    print(data)  # 打印结果


