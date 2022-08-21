# from pymatsci.model import Graphene
# gra = Graphene(15, 0)
# gra.write_vasp('./POCSAR')
# from pymatsci.model import Nanotube
# model = Nanotube(10, 5, 1.42, ['C'], 5)
# model.write_vasp('./POSCAR')

from  pymatsci.model import MagicGraphene
magic = MagicGraphene(10, 5)
magic.write_vasp('./POSCAR')

# from pymatsci.correction import ThermalCorrection
# t = ThermalCorrection(298.15, 101325, True, 3)
# t.free_gas_correction()
# t.printout()


