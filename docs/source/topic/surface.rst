pymatgen创建不同终端的表面

.. code:: python

   from pymatgen.core.surface import generate_all_slabs
   from pymatgen.core.structure import Structure
   from pymatgen.io.vasp import Poscar
   structure = Structure.from_file('CONTCAR')   # 读取文件
   slabs = generate_all_slabs(structure, 2, 1, 0.1, max_normal_search=2)  # 输入结构，最大miller_index，min_slab_size ，min_vacuum_size
   index = 0
   for slab in slabs:
       print(slab.miller_index)
       poscar = Poscar(slab, comment=str(slab.miller_index))
       poscar.write_file(str(index)+str(slab.miller_index)+'.vasp', direct=False)
       index += 1



pymatgen与ASE结构转换

.. code:: python

   from pymatgen.io.ase import AseAtomsAdaptor
   structure = Structure.from_file('CONTCAR')
   atoms = adaptor.get_atoms(slab)
   structure = read('CONTCAR')

ASE创建不对称表面

 .. code:: python

   structure = read('CONTCAR')
   slab = surface(structure, [1, 1, 1], 1, vacuum=2)
   vaccum = 15
   slab.cell[2, 2] = slab.cell[2, 2] + vaccum - 2*2
   # print(slab.get_positions())
   coords = slab.get_positions()
   z = np.around(coords[:, 2], 5)
   unique_z = np.unique(z, axis=0)
   tags_list = []
   index = 1
   for j in range(len(z)):
       for i in range(len(unique_z)):
           if z[j] == unique_z[i]:
               tags_list.append(i+1)
   slab.set_tags(tags_list)
   layers_num = len(unique_z)
   print(layers_num)
   write('POSCAR', slab)