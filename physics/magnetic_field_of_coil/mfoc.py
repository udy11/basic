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

    atol = 1.0e-12
    rtol = 1.0e-12

    mu = 4.0e-7*np.pi
    cr = 0.417
    ci = 1.38e6
    ax = np.array([0.0, 0.0, 1.0])
    xc = np.array([0.0, 0.2, 0.1])
    xr = np.array([0.3, -0.4, 2.0])

    print('----- CHECK 01: TRANSLATION -----')
    print('Checking if output is same if xc and xr are moved by same displacement vector')
    dx = np.random.uniform(-100, 100, size=(3,))
    bm0 = mfoc(mu, cr, ci, ax, xc, xr)
    bm1 = mfoc(mu, cr, ci, ax, xc+dx, xr+dx)
    print(f'B0: {bm0[0]:.15e}, {bm0[1]:.15e}, {bm0[2]:.15e}')
    print(f'After moving by: {dx[0]:.15f}, {dx[1]:.15f}, {dx[2]:.15f}')
    print(f'B1: {bm1[0]:.15e}, {bm1[1]:.15e}, {bm1[2]:.15e}')
    print(f'Passed? {np.isclose(bm0, bm1, atol=atol, rtol=rtol)}')

    print()
    print('----- CHECK 02: ROTATION -----')
    print('Checking if output is same if ax and xr are rotated by same angles')
    a, b, g = np.random.uniform(0.0, 2.0*np.pi, size=(3,))
    ca, cb, cg = np.cos(a), np.cos(b), np.cos(g)
    sa, sb, sg = np.sin(a), np.sin(b), np.sin(g)
    rot0 = np.array([[ca*cb, ca*sb*sg-sa*cg, ca*sb*cg+sa*sg],
                     [sa*cb, sa*sb*sg+ca*cg, sa*sb*cg-ca*sg],
                     [-sb,   cb*sg,          cb*cg]])
    rot1 = np.array([[cb*cg, sa*sb*cg-ca*sg, ca*sb*cg+sa*sg],
                    [cb*sg, sa*sb*sg+ca*cg, ca*sb*sg-sa*cg],
                    [-sb,   sa*cb,          ca*cb]])
    rot = rot0    # choose whether to use rot0 or rot1
    bm0 = mfoc(mu, cr, ci, ax, xc, xr)
    bm1_unrot = mfoc(mu, cr, ci, np.dot(rot, ax), xc, np.dot(rot, xr-xc)+xc)
    bm1 = np.dot(np.linalg.inv(rot), bm1_unrot)    # rotating back the new B vector
    print(f'B0: {bm0[0]:.15e}, {bm0[1]:.15e}, {bm0[2]:.15e}')
    print(f'After rotating by Euler angles: {a:.15f}, {b:.15f}, {g:.15f}')
    print(f'B1: {bm1[0]:.15e}, {bm1[1]:.15e}, {bm1[2]:.15e}')
    print(f'Passed? {np.isclose(bm0, bm1, atol=atol, rtol=rtol)}')

    print()
    print('----- CHECK 03: REFLECTION -----')
    print('Checking if output is same if ax and xr are rotated by same angles')
    rx, ry, rz = [np.ones(3)] * 3
    rx[0], ry[1], rz[2] = -1.0, -1.0, -1.0
    bm0 = mfoc(mu, cr, ci, ax, xc, xr)
    bm1_x = mfoc(mu, cr, ci, ax, xc, (xr-xc)*rx+xc)
    bm1_y = mfoc(mu, cr, ci, ax, xc, (xr-xc)*ry+xc)
    bm1_z = mfoc(mu, cr, ci, ax, xc, (xr-xc)*rz+xc)
    print(f'B0:   {bm0[0]:.15e}, {bm0[1]:.15e}, {bm0[2]:.15e}')
    print(f'After reflecting xr in x, y and z:')
    print(f'B1_x: {bm1_x[0]:.15e}, {bm1_x[1]:.15e}, {bm1_x[2]:.15e}')
    print(f'Passed? {np.isclose(bm0, bm1_x, atol=atol, rtol=rtol)}')
    print(f'B1_y: {bm1_y[0]:.15e}, {bm1_y[1]:.15e}, {bm1_y[2]:.15e}')
    print(f'Passed? {np.isclose(bm0, bm1_y, atol=atol, rtol=rtol)}')
    print(f'B1_z: {bm1_z[0]:.15e}, {bm1_z[1]:.15e}, {bm1_z[2]:.15e}')
    print(f'Passed? {np.isclose(bm0, bm1_z, atol=atol, rtol=rtol)}')
