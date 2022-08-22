
import numpy as np
from pymatgen.core import Structure


class Surface():
    """输出表面吸附结构参数"""
    def __init__(self, index_adorp_atom, num_layer, relaxed_slab):
        structure = Structure.from_file("./CONTCAR")
        coords = structure.cart_coords
        self.index_slab_atom = [i for i in range(coords.shape[0])]
        for i in index_adorp_atom:
            self.index_slab_atom = self.index_slab_atom.remove(i)
        self.num_adorp_atom = len(index_adorp_atom)              
        self.num_layer = num_layer
        self.relaxed_slab = relaxed_slab

        """
        argv:
            num_adorp_atom： 吸附原子数
            num_layer：slab层数
            relaxed_slab：弛豫的表面层数
        
        
        """
    def get_structure_info(self):

        structure = Structure.from_file("./CONTCAR")
        # num_adorp_atom = int(input("请输入吸附原子数："))
        # num_layer = int(input("请输入slab层数："))
        # relaxed_slab = int(input("请输入弛豫的表面层数："))
        coords = structure.cart_coords
        total_atoms = coords.shape[0]  # 总原子数
        slab_atoms = total_atoms - self.num_adorp_atom  #平板原子数
        coords_slab = coords[0:slab_atoms]   # 平板坐标
        coords_adorp_atoms = coords[slab_atoms:, :]  # 吸附原子坐标
        atoms_per_layper = int((total_atoms - num_adorp_atom) / num_layer)  # 每层slab的原子数
        coords_slab = coords_slab[np.argsort(-coords_slab[:, 2])]  # 按z值从大到小排列

        # 求每层原子的平均z值
        z_layer_list = []
        for i in range(num_layer):
            z_layer = 0
            for j in range(atoms_per_layper):
                z_layer += coords_slab[i*atoms_per_layper+j, 2]
            z_layer = z_layer / atoms_per_layper
            print("第{}层原子的z为：{}".format(i+1, z_layer))
            z_layer_list.append(z_layer)
        # 求层间距
        for i in range(relaxed_slab+1):
            dis = z_layer_list[i] - z_layer_list[i+1]
            print("第{}层与第{}层的间距为：{}".format(i+1, i+2, dis))
        # 求吸附原子与表面的距离
        if num_adorp_atom > 0:
            for i in range(num_adorp_atom):
                dis = coords_adorp_atoms[i, 2] - z_layer_list[0]
                print("吸附原子{}与表面的距离为：{}".format(i+1, dis))
        # 求吸附原子之间的距离
        if num_adorp_atom > 1:
            dis = structure.get_distance(-1, -2)
            print("吸附原子间的距离为：%f" % dis)
