# File : setup.py 
  
from distutils.core import setup, Extension 
#name of module 
name  = "cec2020BoundConstrained"
  
#version of module 
version = "1.0"
  
# specify the name of the extension and source files 
# required to compile this 
ext_modules = Extension(name='_cec2020BoundConstrained',sources=["_cec2020BoundConstrained_module.cc", "src/cec20_test_func.cpp", "src/main.cpp"]) 
  
setup(name=name, 
      version=version, 
      ext_modules=[ext_modules])


"""
***COMPILING***
swig -python -c++ -o _cec2020BoundConstrained_module.cc cec2020BoundConstrained.i
python3 setupCec2020BoundConstrained.py build_ext --inplace

***IMPORTING***
import sys
sys.path.insert(0,'nameOfFolder')
import nameOfFile

"""