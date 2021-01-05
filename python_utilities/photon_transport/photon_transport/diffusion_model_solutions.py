"""
Various diffusion model solutions for different number of layers and source functions.
"""

import numpy as np

A_default = 0.14386

def ReflIsoL2(mua1, mua2, mus1r, mus2r, d1, A=A_default):
    """
    Calculate diffuse reflectance from two-layered skin model, the bottom layer
    being semi-infinite.  From "Tissue Parameters Determining the Visual
    Appearance of Normal Skin and Port-Wine Stains", L. O. Svaasand, L. T.
    Norvang et al., Lasers in Medical Science 10, 1995 (complete analytical
    solution to two-layered case in appendix)

    NOTE: Due to numerical inaccuracies, this will fail if mutr**2 * del**2 is
    close to or equal to 1 (i.e. del = 3*D, or mua = musr/2), as the source
    term will blow up. Reproduced by setting mua = 500, musr = 1000.

    Parameters
    ----------
    mua1:
        Absorption in first layer
    mua2:
        Absorption in second layer
    mus1r:
        Reduced scattering coefficient in first layer
    mus2r:
        Reduced scattering coefficient in second layer
    d1:
        Thickness of first layer
    Returns
    -------
    gamma:
        Diffuse reflectance
    """
    D1 = 1/(3*(mus1r + mua1)) #diffusion constant of first layer
    D2 = 1/(3*(mus2r + mua2)) #diffusion constant of second layer

    del1 = np.sqrt(D1/mua1) #optical penetration depth of first layer
    del2 = np.sqrt(D2/mua2) #optical penetration depth of second layer

    #analytical solution to the diffusion equations with continuous boundary conditions at each layer
    fact1 = del1*mus1r*A
    delledd1 = (del1**2*del2/3-del1**2*D2)*np.cosh(d1/del1)
    delledd2 = (del1**3*D2/(3*D1) - del1*del2*D1)*np.sinh(d1/del1)
    ledd1 = (delledd1+delledd2)*(1+del2/(3*D2))
    ledd2 = mus2r/mus1r*del2**2*D1*(del1**2/(9*D1**2)-1)
    ledd3 = del1**2*(D2-del2**2/(9*D2))
    fact2 = ledd1 + (ledd2 + ledd3)*np.exp(-d1/(3*D1))
    fact3 = (del1**2/(9*D1**2) - 1)*(del2/(3*D2)+1)*(D1*del1*(D2+del2*A)*np.cosh(d1/del1)+(D1**2*del2 + D2*del1**2*A)*np.sinh(d1/del1))
    gamma = fact1*fact2*1/fact3
    return gamma

def ReflIsoL3(mua1, mua2, mua3, musr1, musr2, musr3, d1, d2, A=A_default):
    """
    Diffuse reflectance from a three-layered skin model, semi-infinite bottom layer.

    Parameters
    ----------
    mua1:
        Absorption in first layer
    mua2:
        Absorption in second layer
    mua3:
        Absorption in third layer
    musr1:
        Reduced scattering in first layer
    musr2:
        Reduced scattering in second layer
    musr3:
        Reduced scattering in third layer
    d1:
        z-position of the top of the second layer (thickness of first layer)
    d2:
        z-position of the top of the third layer (thickness of first layer + thickness of second layer)
    Returns
    -------
    gamma:
        Diffuse reflectance
    """

    #diffusion constants
    D1 = 1.0/(3.0*(musr1 + mua1));
    D2 = 1.0/(3.0*(musr2 + mua2));
    D3 = 1.0/(3.0*(musr3 + mua3));

    #transport coefficients
    mutr1 = musr1 + mua1;
    u1 = mutr1;
    mutr2 = musr2 + mua2;
    u2 = mutr2;
    mutr3 = musr3 + mua3;
    u3 = mutr3;

    #optical penetration depth
    del1 = np.sqrt(D1/(mua1));
    a1 = del1;
    del2 = np.sqrt(D2/(mua2));
    a2 = del2;
    del3 = np.sqrt(D3/(mua3));
    a3 = del3;

    K1 = del1*del1*musr1/(D1*(1-mutr1*mutr1*del1*del1));
    K2 = del2*del2*musr2/(D2*(1-mutr2*mutr2*del2*del2));
    K3 = del3*del3*musr3/(D3*(1-mutr3*mutr3*del3*del3));
    A2 = -(-2*D2*np.exp(-u1*d1)*np.exp(u2*(-d2+d1))*(K3-K2)*D1*D3*a2-D1*np.exp(-(u1*d1*a2-d1+d2)/a2)*D2*a2*D3*K2+D1*np.exp(-(u1*d1*a2-d1+d2)/a2)*D2*a2*D3*K1+D1*np.exp(-(-d2+d1+u1*d1*a2)/a2)*D2*a2*D3*K1-D1*np.exp(-(-d2+d1+u1*d1*a2)/a2)*D2*a2*D3*K2-np.exp((-d2+d1)/a2)*D3*a2*a2*D1*np.exp(-u1*d1)*(D2*u2*K2-D1*u1*K1)+np.exp(-(-d2+d1)/a2)*D3*a2*a2*D1*np.exp(-u1*d1)*(D2*u2*K2-D1*u1*K1)+D1*np.exp(-(u1*d1*a2-d1+d2)/a2)*D2*D2*a3*K2-D1*np.exp(-(u1*d1*a2-d1+d2)/a2)*D2*D2*a3*K1-D1*np.exp(-(-d2+d1+u1*d1*a2)/a2)*D2*D2*a3*K2+D1*np.exp(-(-d2+d1+u1*d1*a2)/a2)*D2*D2*a3*K1-2*A*a1*D2*np.exp(-u1*d1)*np.exp(u2*(-d2+d1))*(K3-K2)*D3*a2-2*np.exp(-u1*d1-u2*d2+u2*d1)*a2*D2*D2*u2*K2*a3*D1+A*a1*np.exp(-(-d2+d1+u1*d1*a2)/a2)*D2*a2*D3*K1-A*a1*np.exp(-(-d2+d1+u1*d1*a2)/a2)*D2*a2*D3*K2+A*a1*np.exp(-(u1*d1*a2-d1+d2)/a2)*D2*a2*D3*K1-A*a1*np.exp(-(u1*d1*a2-d1+d2)/a2)*D2*a2*D3*K2+A*a1*np.exp(-(-d2+d1)/a2)*a3*D2*a2*np.exp(-u1*d1)*(D2*u2*K2-D1*u1*K1)+A*a1*np.exp((-d2+d1)/a2)*D2*a3*a2*np.exp(-u1*d1)*(D2*u2*K2-D1*u1*K1)+2*np.exp(-u1*d1-u2*d2+u2*d1)*a2*D3*K3*u3*a3*D2*D1+A*a1*np.exp(-(u1*d1*a2-d1+d2)/a2)*D2*D2*a3*K2-A*a1*np.exp(-(u1*d1*a2-d1+d2)/a2)*D2*D2*a3*K1+A*a1*np.exp(-(-d2+d1+u1*d1*a2)/a2)*D2*D2*a3*K1-A*a1*np.exp(-(-d2+d1+u1*d1*a2)/a2)*D2*D2*a3*K2+np.exp((-d2+d1)/a2)*D2*a3*a2*D1*np.exp(-u1*d1)*(D2*u2*K2-D1*u1*K1)+np.exp(-(-d2+d1)/a2)*a3*D2*a2*D1*np.exp(-u1*d1)*(D2*u2*K2-D1*u1*K1)-A*a1*np.exp((-d2+d1)/a2)*D3*a2*a2*np.exp(-u1*d1)*(D2*u2*K2-D1*u1*K1)+A*a1*np.exp(-(-d2+d1)/a2)*D3*a2*a2*np.exp(-u1*d1)*(D2*u2*K2-D1*u1*K1)+2*A*a1*np.exp(-u1*d1-u2*d2+u2*d1)*a2*D3*K3*u3*a3*D2-2*A*a1*np.exp(-u1*d1-u2*d2+u2*d1)*a2*D2*D2*u2*K2*a3-A*a1*a2*D2*np.exp((-d1*a2+d1*a1-d2*a1)/a1/a2)*K1*D3-A*a1*np.exp(-(d1*a2-d2*a1+d1*a1)/a1/a2)*a2*D2*K1*D3+D1*D1*a2*D2*np.exp((-d1*a2+d1*a1-d2*a1)/a1/a2)*K1*u1*a3+D1*a1*D2*D2*np.exp((-d1*a2+d1*a1-d2*a1)/a1/a2)*K1*u1*a3+A*D1*a2*D2*np.exp((-d1*a2+d1*a1-d2*a1)/a1/a2)*K1*a3+D1*D1*np.exp(-(d1*a2-d2*a1+d1*a1)/a1/a2)*a2*D2*K1*u1*a3+A*D1*np.exp(-(d1*a2-d2*a1+d1*a1)/a1/a2)*a2*D2*K1*a3-D1*a1*np.exp(-(d1*a2-d2*a1+d1*a1)/a1/a2)*D2*D2*K1*u1*a3+A*D1*np.exp(-(d1*a2-d2*a1+d1*a1)/a1/a2)*a2*a2*K1*D3+D1*D1*np.exp(-(d1*a2-d2*a1+d1*a1)/a1/a2)*a2*a2*K1*u1*D3-D1*D1*a2*a2*np.exp((-d1*a2+d1*a1-d2*a1)/a1/a2)*K1*u1*D3-A*D1*a2*a2*np.exp((-d1*a2+d1*a1-d2*a1)/a1/a2)*K1*D3+A*a1*D2*D2*np.exp((-d1*a2+d1*a1-d2*a1)/a1/a2)*K1*a3-A*a1*np.exp(-(d1*a2-d2*a1+d1*a1)/a1/a2)*D2*D2*K1*a3-D1*a1*a2*D2*np.exp((-d1*a2+d1*a1-d2*a1)/a1/a2)*K1*u1*D3-D1*a1*np.exp(-(d1*a2-d2*a1+d1*a1)/a1/a2)*a2*D2*K1*u1*D3)*a1/(np.exp(-(-d1*a2+d1*a1-d2*a1)/a1/a2)*D1*a1*D2*D2*a3+np.exp(-(-d1*a2+d1*a1-d2*a1)/a1/a2)*D1*D1*a2*D2*a3+np.exp(-(-d1*a2+d1*a1-d2*a1)/a1/a2)*A*a1*a1*D2*D2*a3+np.exp(-(d1*a2-d2*a1+d1*a1)/a1/a2)*D1*a1*D2*D2*a3-np.exp(-(d1*a2-d2*a1+d1*a1)/a1/a2)*D1*D1*a2*D2*a3-np.exp(-(d1*a2-d2*a1+d1*a1)/a1/a2)*A*a1*a1*D2*D2*a3+np.exp((-d1*a2+d1*a1-d2*a1)/a1/a2)*A*a1*a1*D2*D2*a3-np.exp((-d1*a2+d1*a1-d2*a1)/a1/a2)*D1*D1*a2*D2*a3-np.exp((-d1*a2+d1*a1-d2*a1)/a1/a2)*D1*a1*D2*D2*a3-np.exp((d1*a2-d2*a1+d1*a1)/a1/a2)*A*a1*a1*D2*D2*a3+np.exp((d1*a2-d2*a1+d1*a1)/a1/a2)*D1*D1*a2*D2*a3-np.exp((d1*a2-d2*a1+d1*a1)/a1/a2)*D1*a1*D2*D2*a3+np.exp(-(-d1*a2+d1*a1-d2*a1)/a1/a2)*A*a1*a1*D2*a2*D3+np.exp(-(-d1*a2+d1*a1-d2*a1)/a1/a2)*D1*a1*D2*a2*D3+np.exp(-(-d1*a2+d1*a1-d2*a1)/a1/a2)*A*a1*D1*a2*a2*D3+np.exp((d1*a2-d2*a1+d1*a1)/a1/a2)*D1*a1*D2*a2*D3+np.exp((d1*a2-d2*a1+d1*a1)/a1/a2)*A*a1*a1*D2*a2*D3-np.exp((-d1*a2+d1*a1-d2*a1)/a1/a2)*A*a1*a1*D2*a2*D3+np.exp(-(d1*a2-d2*a1+d1*a1)/a1/a2)*A*a1*D1*a2*a2*D3-np.exp(-(d1*a2-d2*a1+d1*a1)/a1/a2)*A*a1*a1*D2*a2*D3-np.exp((-d1*a2+d1*a1-d2*a1)/a1/a2)*A*a1*D1*a2*a2*D3+np.exp(-(d1*a2-d2*a1+d1*a1)/a1/a2)*D1*a1*D2*a2*D3-np.exp((d1*a2-d2*a1+d1*a1)/a1/a2)*A*a1*D1*a2*a2*D3+np.exp((-d1*a2+d1*a1-d2*a1)/a1/a2)*D1*a1*D2*a2*D3+np.exp((-d1*a2+d1*a1-d2*a1)/a1/a2)*A*a1*D1*a2*D2*a3+np.exp((d1*a2-d2*a1+d1*a1)/a1/a2)*A*a1*D1*a2*D2*a3+np.exp(-(d1*a2-d2*a1+d1*a1)/a1/a2)*A*a1*D1*a2*D2*a3+np.exp(-(-d1*a2+d1*a1-d2*a1)/a1/a2)*A*a1*D1*a2*D2*a3+np.exp(-(-d1*a2+d1*a1-d2*a1)/a1/a2)*D1*D1*a2*a2*D3-np.exp((d1*a2-d2*a1+d1*a1)/a1/a2)*D1*D1*a2*a2*D3-np.exp(-(d1*a2-d2*a1+d1*a1)/a1/a2)*D1*D1*a2*a2*D3+np.exp((-d1*a2+d1*a1-d2*a1)/a1/a2)*D1*D1*a2*a2*D3);
    A1 = (A2*(D1/del1-A)-K1*(A+D1*mutr1))/(D1/del1+A);
    j = K1*D1*mutr1 + D1/del1*A1 - A2*D1/del1;
    r = -j;
    return r


def ReflIsoL1(mua, musr, A=A_default):
    """
    Calculate reflectance from a one-layered, semi-infinite skin model.

    From Bjorgan 2013, master thesis, p. 23
    """
    D = 1/(3*(musr + mua)) #diffusion constant
    pen_depth = np.sqrt(D/mua) #optical penetration depth

    return musr*A*pen_depth*pen_depth/((pen_depth/(3*D) + 1.0)*(D + pen_depth*A))

def ReflE2L2(mua1, mua2, mus1, mus2, g1, g2, d1, A=A_default):
    """
    Two layers, delta-eddington.
    """
    #reduced transport coefficients
    mutrr1 = mus1*(1-g1**2)+mua1
    mutrr2 = mus2*(1-g2**2)+mua2

    mutr1 = mus1*(1-g1)+mua1
    mutr2 = mus2*(1-g2)+mua2

    #diffuse transport coefficients
    D1 = 1/(3*(mus1*(1-g1)+mua1))
    D2 = 1/(3*(mus2*(1-g2)+mua2))

    #optical penetration depths
    del1 = np.sqrt(D1/mua1)
    del2 = np.sqrt(D2/mua2)

    #source function factors
    K1 = mus1*(1-g1**2)*(1/D1+3*mutrr1*g1/(1+g1))*del1**2/(1-mutrr1**2*del1**2)
    K2 = mus2*(1-g2**2)*(1/D2+3*mutrr2*g2/(1+g2))*del2**2/(1-mutrr2**2*del2**2)*np.exp(-mutrr1*d1)

    #diffuse reflectance
    gamma = -A * (((D1 * K1 * mutrr1 * mutr1 - mus1 * (-1 + g1) * g1) * del1 + mutr1 * D1 * K1) * mutr2 * (D1 * del2 - D2 * del1) * np.exp(-(d1 / del1)) + mutr2 * ((D1 * K1 * mutrr1 * mutr1 - mus1 * (-1 + g1) * g1) * del1 - mutr1 * D1 * K1) * (D1 * del2 + D2 * del1) * np.exp((d1 / del1)) - 0.2e1 * del1 * D1 * (mutr2 * (D1 * mutr1 * del2 * K1 * mutrr1 - g1 * mus1 * (-1 + g1) * del2 - D2 * mutr1 * K1) * np.exp(-(mutrr1 * d1)) - mutr1 * (-g2 * mus2 * del2 * (-1 + g2) * np.exp(-(mutrr2 * d1)) + K2 * D2 * mutr2 * (-1 + mutrr2 * del2)))) / mutr1 / mutr2 / ((D1 * del2 - D2 * del1) * (-D1 + A * del1) * np.exp(-(d1 / del1)) + (D1 * del2 + D2 * del1) * np.exp((d1 / del1)) * (D1 + A * del1))
    return gamma

def ReflE2L3(mua1, mua2, mua3, mus1, mus2, mus3, g1, g2, g3, d1, d2, A=A_default):
    """
    Three layers, delta-eddington.
    """
    #reduced transport coefficients
    mutrr1 = mus1*(1-g1**2)+mua1;
    mutrr2 = mus2*(1-g2**2)+mua2;
    mutrr3 = mus3*(1-g3**2)+mua3;

    mutr1 = mus1*(1-g1)+mua1;
    mutr2 = mus2*(1-g2)+mua2;
    mutr3 = mus3*(1-g3)+mua3;

    #diffuse transport coefficients
    D1 = 1/(3*(mus1*(1-g1)+mua1));
    D2 = 1/(3*(mus2*(1-g2)+mua2));
    D3 = 1/(3*(mus3*(1-g3)+mua3));

    #optical penetration depths
    del1 = np.sqrt(D1/mua1);
    del2 = np.sqrt(D2/mua2);
    del3 = np.sqrt(D3/mua3);

    #factors in the source functions
    K1 = mus1*(1-g1**2)*(1/D1+3*mutrr1*g1/(1+g1))*del1**2/(1-mutrr1**2*del1**2);
    K2 = mus2*(1-g2**2)*(1/D2+3*mutrr2*g2/(1+g2))*del2**2/(1-mutrr2**2*del2**2)*np.exp(-mutrr1*d1);
    K3 = mus3*(1-g3**2)*(1/D3+3*mutrr3*g3/(1+g3))*del3**2/(1-mutrr3**2*del3**2)*np.exp(-mutrr1*d1)*np.exp(-mutrr2*(d2-d1));

    #diffuse reflectance
    gamma = -((np.exp(d2 / del2) * (D2 * del3 + D3 * del2) * (-D2 * del1 + D1 * del2) * ((D1 * mutr1 * K1 * mutrr1 - (mus1 * (-1 + g1) * g1)) * del1 + D1 * mutr1 * K1) * mutr2 * mutr3 * np.exp(-d1 / del1) - 0.4e1 * D1 * D2 * np.exp(d1 / del2) * del2 * mutr1 * del1 * mutr2 * mutr3 * K2 * (D2 * del3 * mutrr2 - D3) * np.exp(mutrr2 * (-d2 + d1)) + 0.4e1 * D2 * (g2 * del3 * mus2 * mutr3 * (-1 + g2) * np.exp(-mutrr2 * d2) + (-g3 * del3 * mus3 * (-1 + g3) * np.exp(-mutrr3 * d2) + D3 * K3 * mutr3 * (del3 * mutrr3 - 0.1e1)) * mutr2) * D1 * del1 * del2 * mutr1 * np.exp(d1 / del2) + np.exp(d2 / del2) * (D2 * del3 + D3 * del2) * ((D1 * del2 + D2 * del1) * ((D1 * mutr1 * K1 * mutrr1 - (mus1 * (-1 + g1) * g1)) * del1 - D1 * mutr1 * K1) * mutr2 * np.exp(d1 / del1) - 0.2e1 * ((D1 * mutr1 * del2 * K1 * mutrr1 - g1 * mus1 * (-1 + g1) * del2 - K1 * mutr1 * D2) * mutr2 * np.exp(-mutrr1 * d1) - (-g2 * mus2 * del2 * (-1 + g2) * np.exp(-mutrr2 * d1) + D2 * K2 * mutr2 * (mutrr2 * del2 - 0.1e1)) * mutr1) * D1 * del1) * mutr3) * np.exp(-d1 / del2) + np.exp(-d2 / del2) * ((D1 * del2 + D2 * del1) * ((D1 * mutr1 * K1 * mutrr1 - (mus1 * (-1 + g1) * g1)) * del1 + D1 * mutr1 * K1) * mutr2 * np.exp(-d1 / del1) + ((D1 * mutr1 * K1 * mutrr1 - (mus1 * (-1 + g1) * g1)) * del1 - D1 * mutr1 * K1) * (-D2 * del1 + D1 * del2) * mutr2 * np.exp(d1 / del1) - 0.2e1 * D1 * del1 * ((D1 * mutr1 * del2 * K1 * mutrr1 - g1 * mus1 * (-1 + g1) * del2 + K1 * mutr1 * D2) * mutr2 * np.exp(-mutrr1 * d1) - (-g2 * mus2 * del2 * (-1 + g2) * np.exp(-mutrr2 * d1) + D2 * K2 * mutr2 * (mutrr2 * del2 + 0.1e1)) * mutr1)) * np.exp(d1 / del2) * mutr3 * (D2 * del3 - D3 * del2)) * A / (((-D2 * del1 + D1 * del2) * (-D1 + A * del1) * np.exp(-d1 / del1) + (D1 * del2 + D2 * del1) * np.exp(d1 / del1) * (A * del1 + D1)) * np.exp(d2 / del2) * (D2 * del3 + D3 * del2) * np.exp(-d1 / del2) + np.exp(-d2 / del2) * ((D1 * del2 + D2 * del1) * (-D1 + A * del1) * np.exp(-d1 / del1) + (-D2 * del1 + D1 * del2) * np.exp(d1 / del1) * (A * del1 + D1)) * np.exp(d1 / del2) * (D2 * del3 - D3 * del2)) / mutr2 / mutr1 / mutr3;
    return gamma


