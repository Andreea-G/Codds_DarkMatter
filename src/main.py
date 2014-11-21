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

options = {0 : zero,
           1 : sqr,
           4 : sqr,
           9 : sqr,
           2 : lambda x: x**2,
           3 : prime,
           5 : prime,
           7 : prime,
           'x' : lambda : "ha",
}
 
        
class Experiment:
    def __init__(self, expername):
        module = import_file(INPUT_DIR + expername + ".py")
        self.name = expername
        self.var = module.myvar
        self.fnc = options[self.var]
    def getvar(self):
        return self.var
    def myfnc(self):
        print options[self.var]()
        
        
def main():
    exper = "LUX"
    cl1 = Experiment(exper)
    print 'var = ', cl1.var
    print 'name = ', cl1.name
    print 'getvar = ', cl1.getvar()
    print 'myfunc: '
    cl1.myfnc()
    print cl1.fnc()    

    
if __name__ == '__main__':
    main()

