"""
Various derivatives of reflectance model in photon_transport.forward_model,
used for simple_inverse_model.
"""

import numpy as np
from photon_transport import forward_model

def ReflIsoL2DerivMuad(de, mua1, mua2, musr1, musr2):
    """
    Calculate the derivative of the simulated two-layered reflectance with
    respect to the absorption in the lower layer.

    Parameters
    ----------
    de: float
        Thickness of upper layer
    mua1: float
        Absorption coefficient in upper layer
    mua2: float
        Absorption coefficient in lower layer
    musr1: float
        Reduced scattering coefficient in upper layer
    musr2: float
        Reduced scattering coefficient in lower layer
    Returns
    -------
    derivative: float
        Derivative
    """
    D1 = 1.0/(3.0*(musr1 + mua1))
    musr2dmusr1 = musr2/musr1
    A = forward_model.A
    div13 = 1.0/3.0

    #fixed parameters related to epidermis
    del1 = np.sqrt(np.divide(D1,mua1))
    sinhval = np.sinh(np.divide(de, del1))
    coshval = np.cosh(np.divide(de, del1))
    expval = np.exp(-np.divide(de,D1)*div13)

    #dermis parameters
    D2 = np.divide(1.0,3.0*(musr2 + mua2))
    del2 = np.sqrt(np.divide(D2,mua2))

    #calculate reflectance directly
    f1 = (del1*del1*del2*div13-del1*del1*D2)*coshval+(del1*del1*del1*np.divide(D2,D1)*div13 - del1*del2*D1)*sinhval
    f2 = 1.0+np.divide(del2, D2)*div13
    f3 = (musr2dmusr1*del2*del2*D1*(del1*np.divide(del1,D1*D1)*div13*div13-1.0)+del1*del1*(D2-del2*np.divide(del2,D2)*div13*div13))*expval
    f4 = del1*np.divide(del1,D1*D1)*div13*div13 - 1.0
    f5 = np.divide(del2,D2)*div13+1.0
    f6 = D1*del1*(D2+del2*A)*coshval+(D1*D1*del2 + D2*del1*del1*A)*sinhval
    num = del1*musr1*A*(f1*f2+f3)
    denom = f4*f5*f6
    reflectance = np.divide(num, denom)

    #calculate the derivative with respect to muad
    dD2dmuad = -3.0*D2*D2
    ddel2dmuad = (dD2dmuad*mua2-D2)*np.divide(1.0, mua2*mua2)*np.divide(1.0, 2.0*del2)
    df2dmuad = div13*(ddel2dmuad*D2 - del2*dD2dmuad)*np.divide(1.0, D2*D2)
    derivative = (del1*musr1*A*((coshval*(del1*del1*div13*ddel2dmuad - del1*del1*dD2dmuad) + sinhval*(del1*del1*del1*np.divide(1.0, D1)*div13*dD2dmuad - del1*D1*ddel2dmuad))*f2 + f1*df2dmuad + expval*(musr2dmusr1*2.0*del2*ddel2dmuad*D1*(np.divide(del1*del1,9.0*D1*D1)-1.0) + del1*del1*(dD2dmuad + np.divide(1.0, 9.0*mua2*mua2))))*denom - (f4*(df2dmuad*f6 + D1*del1*(dD2dmuad + ddel2dmuad*A)*coshval + (D1*D1*ddel2dmuad + dD2dmuad*del1*del1*A)*sinhval*f5)*num))*np.divide(1.0, denom*denom)
    return reflectance, derivative

def ReflIsoL2DerivMuae(d1, mua1, mua2, musr1, musr2):
    """
    Calculate derivative of simulated two-layered reflectance with respect to the absorption in the upper layer.
    See ReflIsoL2DerivMuad.
    """
    A = forward_model.A

    #fixed parameters related to dermis
    D2 = np.divide(1.0,3.0*(musr2 + mua2))
    del2 = np.sqrt(np.divide(D2,mua2))

    f2 = 1.0 + np.divide(del2 ,(D2 * 3.0))
    f6 = D2-del2*np.divide(del2, (D2*9.0))

    #epidermis parameters
    D1 = np.divide(1.0,3.0*(musr1 + mua1))
    del1 = np.sqrt(np.divide(D1,mua1))
    div1D1D1 = np.divide(1.0, D1*D1)

    coshval = np.cosh(np.divide(d1, del1))
    sinhval = np.sinh(np.divide(d1, del1))

    #reflectance
    f1 = (del1*del1 * np.divide(del2 , 3.0) - del1*del1 * D2) * coshval + (np.power(del1, 3.0) * np.divide(D2 , D1*3.0) - del1 * del2 * D1) * sinhval
    f3 = np.divide(musr2, musr1) * del2 * del2 * D1 * (del1*del1 * np.divide(div1D1D1,9.0) - 1.0) + del1*del1 * f6
    f4 = np.exp(-np.divide(d1 , (D1 *3.0)))
    f5 = del1*del1 * np.divide(div1D1D1,9.0) - 1.0
    f7 = D1 * del1 * (D2 + del2 * A) * coshval + (del2 * D1*D1 + D2 * del1*del1 * A) * sinhval
    fact = np.divide((f1 * f2 + f3 * f4) , (f5 * f2 * f7))
    reflectance = del1 * musr1 * A * fact

    #calculate derivative with respect to muae
    dD1d1 = -3.0*D1*D1
    ddel1d1 = (dD1d1*mua1-D1)*np.divide(1.0, mua1*mua1)*np.divide(1.0, 2.0*del1)


    df1d1 = (np.divide(2.0 , 3.0) * del1 * del2 * ddel1d1 - 2.0 * del1 * D2 * ddel1d1) * coshval - (del1*del1 * np.divide(del2 , 3.0) - del1*del1 * D2) * sinhval * d1 * np.divide(1.0, del1*del1) * ddel1d1 + (del1*del1 * np.divide(D2 , D1) * ddel1d1 - np.power(del1, 3.0) * D2 * (div1D1D1) * np.divide(dD1d1 , 3.0) - ddel1d1 * del2 * D1 - del1 * del2 * dD1d1) * sinhval - (np.power(del1, 3.0) * np.divide(D2 , (D1 * 3.0)) - del1 * del2 * D1) * coshval * d1 * np.divide(1.0, del1*del1) * ddel1d1
    df3d1 = np.divide(musr2 , musr1) * del2 * del2 * dD1d1 * (del1*del1 * np.divide(div1D1D1,9.0) - 1.0) + np.divide(musr2 , musr1) * del2 * del2 * D1 * (np.divide(2.0 , 9.0) * del1 * (div1D1D1) * ddel1d1 - np.divide(2.0 , 9.0) * del1*del1 * np.power(D1, -3.0) * dD1d1) + 2.0 * del1 * (f6) * ddel1d1
    df4d1 = d1 * (div1D1D1) * dD1d1 * f4 *np.divide(1.0, 3.0)
    df5d1 = np.divide(2.0 , 9.0) * del1 * (div1D1D1) * ddel1d1 - np.divide(2.0 , 9.0) * del1*del1 * np.power(D1, -3.0) * dD1d1
    df7d1 = dD1d1 * del1 * (D2 + del2 * A) * coshval + D1 * ddel1d1 * (D2 + del2 * A) * coshval - np.divide(D1 , del1) * (D2 + del2 * A) * sinhval * d1 * ddel1d1 + (2.0 * del2 * D1 * dD1d1 + 2.0 * D2 * del1 * A * ddel1d1) * sinhval - (del2 * D1*D1 + D2 * del1*del1 * A) * coshval * d1 * np.divide(1.0, del1*del1) * ddel1d1

    derivative = ddel1d1 * musr1 * A * fact + del1 * musr1 * A * (df1d1 * f2 + df3d1 * f4 + f3 * df4d1) *np.divide(1.0, (f5 * f2 * f7)) - del1 * musr1 * A * fact * np.divide(1.0, f5) * df5d1 - del1 * musr1 * A * fact * np.divide(1.0, f7) * df7d1

    return reflectance, derivative
