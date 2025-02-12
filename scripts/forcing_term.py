import numpy as np

def forcing_term(x, y, z, t, Re):
    f = np.zeros(3)
    f[0] = -sin(t)**2*sin(x)*sin(y)**2*sin(z)**2*cos(x) + sin(t)**2*sin(x)*sin(z)**2*cos(x)*cos(y)**2 + 2*sin(t)**2*sin(x)*cos(x)*cos(y)**2*cos(z)**2 + sin(x)*sin(z)*cos(t)*cos(y) - sin(x)*sin(z)*cos(y) + 3*sin(t)*sin(x)*sin(z)*cos(y)/Re
    f[1] = -sin(t)**2*sin(x)**2*sin(y)*sin(z)**2*cos(y) + sin(t)**2*sin(y)*sin(z)**2*cos(x)**2*cos(y) + 2*sin(t)**2*sin(y)*cos(x)**2*cos(y)*cos(z)**2 + sin(y)*sin(z)*cos(t)*cos(x) - sin(y)*sin(z)*cos(x) + 3*sin(t)*sin(y)*sin(z)*cos(x)/Re
    f[2] = -2*sin(t)**2*sin(x)**2*sin(z)*cos(y)**2*cos(z) - 2*sin(t)**2*sin(y)**2*sin(z)*cos(x)**2*cos(z) - 4*sin(t)**2*sin(z)*cos(x)**2*cos(y)**2*cos(z) + 2*cos(t)*cos(x)*cos(y)*cos(z) + cos(x)*cos(y)*cos(z) + 6*sin(t)*cos(x)*cos(y)*cos(z)/Re
    return f
