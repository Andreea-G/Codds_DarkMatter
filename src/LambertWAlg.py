"""
Makeshift LambertW Algorithm -- B/c standard one doesn't work on Hoffman2 cluster for whatever stupid reason...
"""
import math

eps = 0.00000001 # max error allowed
def w0Lambert(x,k): # Lambert W function using Newton's method
    if k!=-1:
        w=1
        v=100
    elif k==-1:
        w=-2
        v=-100
    

    
    while True:
        v=w
        ew = math.exp(w)
        f=w*ew-x
        wNew = w - f/((ew*(w+1.)-(w+2.)*f/(2.*w+2.)))
        if abs(w - wNew) <= eps: break
        w = wNew
    return w


