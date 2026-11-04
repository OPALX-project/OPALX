"""Exact uniform-B orbit moments for the physical-volume uniform disk."""
import numpy as np
import pandas as pd
from numpy.polynomial.legendre import leggauss
R,b,p0=2.,.02,.001/np.sqrt(1-.001**2)
def evaluate(s,n=256):
    x,w=leggauss(n);theta=(x+1)*np.pi/2
    a=b*np.cos(theta)
    weights=w*np.sin(theta)**2*(1+a/R);weights/=weights.sum()
    rows=[]
    for distance in s:
        phi=distance/R
        radius=np.sqrt(R*R+2*R*a*np.cos(phi)+a*a)
        u=(2*R*a*np.cos(phi)+a*a)/(radius+R)
        pu=-p0*a*np.sin(phi)/radius
        q=u-np.sum(weights*u);p=pu-np.sum(weights*pu)
        v=np.sum(weights*q*q)
        residual=p-np.sum(weights*q*p)/v*q if v else p
        eps=np.sqrt(v*np.sum(weights*residual**2))
        rows.append((distance,eps))
    return pd.DataFrame(rows,columns=['s_m','emit_n_u_m'])
