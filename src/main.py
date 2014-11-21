# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 22:52:11 2014

@author: Andreea
"""
from __future__ import absolute_import

INPUT_DIR = "Data/"
OUTPUT_DIR = "Output/"

def import_file(full_path_to_module):
    import os, sys
    directory, module_name = os.path.split(full_path_to_module)
    module_name = os.path.splitext(module_name)[0]
        
    path = list(sys.path)
    sys.path.insert(0, directory)
    try:
        module = __import__(module_name)
        return module
    finally:
        sys.path[:] = path # restore 
        
        
        
def zero():
    print "You typed zero.\n"
 
def sqr():
    print "n is a perfect square\n"
 
def even():
    print "n is an even number\n"
 
def prime():
    print "n is a prime number\n"

def myfunc2(exp):
    print "Calling myfunc2"
    print exp.name          

options = {0 : zero,
           1 : sqr,
           4 : sqr,
           9 : sqr,
           2 : lambda x: x**2,
           3 : prime,
           5 : prime,
           7 : prime,
           'x' : lambda : "ha",
           'me' : myfunc2,
}
 

def CrossSectionFactors_SI(exper):
    print exper.name, "SI"

def CrossSectionFactors_SD66(exper):
    print exper.name, "SID66"

def CrossSectionFactors_SD44(exper):
    print exper.name, "SD44"

CrossSectionFactors_options = {'SI' : CrossSectionFactors_SI,
           'SD66' : CrossSectionFactors_SD66,
           'SD44' : CrossSectionFactors_SD44,
}

class Experiment:
    def __init__(self, expername, scattering_type):
        module = import_file(INPUT_DIR + expername + ".py")
        self.name = expername
        self.scattering_type = scattering_type
        self.var = module.myvar
        self.fnc = options[self.var]
        self.energy_resolution_type = module.energy_resolution_type
        self.target_nuclide_AZC_list = module.target_nuclide_AZC_list
        self.target_nuclide_JSnSp_list = module.target_nuclide_JSnSp_list

    def getvar(self):
        return self.var

    def myfnc(self):
        print options[self.var](self)
        
    def CrossSectionFactors(self):
        CrossSectionFactors_options[self.scattering_type](self)
              
                
        
def main():
    exper = "superCDMS"
    scattering_type = 'SD66'
    
    
    cl1 = Experiment(exper, scattering_type)
    print 'var = ', cl1.var
    print 'name = ', cl1.name
    print 'getvar = ', cl1.getvar()
    print 'myfunc: '
    cl1.myfnc()   
    print cl1.target_nuclide_AZC_list
    print cl1.target_nuclide_JSnSp_list
    print cl1.energy_resolution_type    
    cl1.CrossSectionFactors()
    
    #myfunc2(cl1)    
    
if __name__ == '__main__':
    main()

