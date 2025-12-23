# Last Updated: 23-Dec-2025
# Udaya Maurya (udaya_cbscients@yahoo.com, telegram: https://t.me/udy11)
# Source: https://github.com/udy11, https://gitlab.com/udy11
# Function to compute magnetic field due to
# current carrying coil at any point in space

# ALL YOU NEED TO DO:
# Call mfoc() with:
#   INPUT:
#     mu = permeability constant of the medium (float)
#     cr = coil radius (float)
#     ci = coil current (float, sign sensitive)
#     ax = direction in which coil axis is pointing (1d-numpy array of size 3, can be unnormalized)
#     xc = co-ordinates of the center of coil (1d-numpy array of size 3)
#     xr = co-ordinates at which magnetic field is to be calculated (1d-numpy array of size 3)
#   OUTPUT:
#     bm = components of magnetic field vector (1d-numpy array of size 3)

# NOTE:
# All numbers and calculations in SI Units
# Please read the attached pdf file for the expression of magnetic field
# A positive ci will be considered anticlockwise with respect to axis direction ax
# ax, xc, xr and bm have 3 elements each representing x, y and z
#   components respectively
# Magnetic field blows up at the coil's position, make sure xr is not
#   very near the coil

import numpy as np
import scipy.special

def mfoc(mu, cr, ci, ax, xc, xr):
    ''' (float, float, float, array, array, array) -> array
        Function to compute magnetic field vector in space due to current carrying coil
    '''
    bm = np.zeros(3)
    x2 = np.zeros(3)
    ax12 = ax[0] * ax[0] + ax[1] * ax[1]
    ax123s = np.sqrt(ax12 + ax[2] * ax[2])
    cth = ax[2] / ax123s
    sth = np.sqrt(ax12) / ax123s
    ph = np.atan2(ax[1],ax[0])
    cph = np.cos(ph)
    sph = np.sin(ph)
    x1 = xr - xc
    x1cs = x1[0] * cph + x1[1] * sph
    x2[0] = x1cs * cth - x1[2] * sth
    x2[1] = -x1[0] * sph + x1[1] * cph
    x2[2] = x1cs * sth + x1[2] * cth
    x2r = np.sqrt(x2[0] * x2[0] + x2[1] * x2[1])
    if x2r < 2.22044604925032e-16:
        tv1 = mu * ci * cr * cr * 0.5 * (1.0 / np.sqrt(cr * cr + x2[2] * x2[2]))**3
        bm[0] = tv1 * cph * sth
        bm[1] = tv1 * sth * sph
        bm[2] = tv1 * cth
    else:
        rt = np.zeros((3, 3))
        x2t = np.atan2(x2[1],x2[0])
        ck2 = 4.0 * cr * x2r / ((x2r + cr)**2 + x2[2]**2)
        if ck2 == 1:
            ek = 1.7976931348623157e+308
        else:
            ek = scipy.special.ellipk(ck2)
        ee = scipy.special.ellipe(ck2)
        tv1 = mu * ci * np.sqrt(ck2) * 0.07957747154594766788 / x2r / np.sqrt(cr * x2r)
        tv2 = 0.5 * ee / (1.0 - ck2)
        tv3 = tv1 * (tv2 * (2.0 - ck2) - ek) * x2[2]
        bm[0] = tv3 * np.cos(x2t)
        bm[1] = tv3 * np.sin(x2t)
        bm[2] = tv1 * (tv2 * (ck2 * (x2r + cr) - 2.0 * x2r) + ek * x2r)
        rt[0,0] = cth * cph
        rt[0,1] = -sph
        rt[0,2] = sth * cph
        rt[1,0] = cth * sph
        rt[1,1] = cph
        rt[1,2] = sth * sph
        rt[2,0] = -sth
        #rt[2,1] = 0.0
        rt[2,2] = cth
        bm = np.dot(rt, bm)
    return bm

if __name__ == '__main__':

    print(mfoc(4.0e-7*np.pi, 1.0, 7.1e5, np.array([0, 0, 1]), np.array([0, 0, 0.1]), np.array([1.39, 0, -0.2])))
