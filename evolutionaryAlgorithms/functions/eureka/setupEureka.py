# File : setup.py 
  
from distutils.core import setup, Extension 
#name of module 
name  = "eureka"
  
#version of module 
version = "1.0"
  
# specify the name of the extension and source files 
# required to compile this 
ext_modules = Extension(name='_eureka',sources=["_eureka_module.cc","src/EurekaOptimaException.cpp","src/Problem.cpp","src/TrussBarStructureStaticProblem.cpp","src/TrussBarStructureStaticSimulator.cpp","src/F101Truss10Bar.cpp","src/F103Truss25Bar.cpp","src/F105Truss60Bar.cpp","src/F107Truss72Bar.cpp","src/F109Truss942Bar.cpp"]) 
  
setup(name=name, 
      version=version, 
      ext_modules=[ext_modules]) 


"""
***COMPILING***
swig -python -c++ -o _eureka_module.cc eureka.i
python3 setupEureka.py build_ext --inplace

***IMPORTING***
import sys
sys.path.insert(0,'nameOfFolder')
import nameOfFile

"""