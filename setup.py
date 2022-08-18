from setuptools import setup, find_packages

setup(name="pymatsci", version='1.0.1', description='Python Materials Science', packages=find_packages(), author='C.L. Qin',
      author_email='clqin@foxmail.com', long_description='Pymatsci (Python Materials Science) is a robust, open-source Python library for materials analysis.', 
      url='https://github.com/stupid-cloud/Pymatsci.git',
      license='MIT License (MIT)', install_requires=['numpy', 'pymatgen'], python_requires='>=3')