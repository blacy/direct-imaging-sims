import numpy as np

def fit_func(w, a, b, c, d):
    s = np.sqrt(1-w)
    return a*(1-b*s)*(1.0-s)/((1. + c*s)*(1 + d*w))

def fit_func_2(w, a, b, c, d):
    s = np.sqrt(1-w)
    return a*(1-b*s)*(1.0-s)/((1. + c*s)*(1+d*s))

def isofit(omega):
    a = 0.6896
    b = -0.4380
    c = 1.2966
    d = 0.7269
    return fit_func_2(omega, a, b, c, d)

def rayfit(omega):
    a = 0.8190
    b = -0.1924
    c = 1.5946
    d = 0.0886
    return fit_func(omega, a, b, c, d)

def vecfit(omega):
    a = 0.6901
    b = -0.4664
    c = 1.7095
    d = -0.1345
    return fit_func(omega, a, b, c, d)
