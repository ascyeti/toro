from sympy import *
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt


def compute_cos(x):
    y = symbols('y')
    f = cos(y)
    return f.evalf(subs={y: x})
def compute_sin(x):
    return compute_cos(x-(pi/2))

def costruisci_toro(r, R, r_density, theta_density):
    x, y = symbols("x y")
    f = sqrt(r**2-(sqrt(x**2+y**2)-R)**2)
    # r è il raggio della circonferenza che ruota attorno all'asse z
    # R è il raggio della circonferenza che passa per tutti i centri delle circonferenze di raggio r

    xyz = [] # array dei punti del toro

    for i in range(r_density):
        for j in range(theta_density):

            r_temp = R - r + i * (2 * r) / r_density
            theta_temp = 2 * pi * j / theta_density

            xp = float(r_temp * cos(theta_temp))
            yp = float(r_temp * sin(theta_temp))
            zp = f.subs([(x, xp), (y, yp)])
            if (not zp.is_real) or zp < 10**(-4): #Messo per bilanciare l'appossimazione
                zp = 0

            xyz.append([xp, yp, zp]) # Poi li sbatto dentro l'array

    print(xyz)

    xyz_interno = []
    xyz_esterno = []
    for p in xyz:
        if p[0]**2 + p[1]**2 < R**2:
            xyz_interno.append(p)
        elif p[0]**2 + p[1]**2 > R**2:
            xyz_esterno.append(p)
        else:
            xyz_interno.append(p)

    xyz_interno = sorted(xyz_interno, key=lambda tup: tup[2])
    matrice_p_interni = np.zeros((int(len(xyz_interno)/theta_density), theta_density, 3))
    matrice_p_esterni = np.zeros((int(len(xyz_esterno) / theta_density), theta_density, 3))

    i = 0
    j = 0

    for p in xyz_interno:
        if xyz_interno.index(p) == 0:
            matrice_p_interni[0][0] = p
            j += 1
        elif j < 20: #p[2] - xyz_interno[xyz_interno.index(p) - 1][2] < 10**(-3):
            matrice_p_interni[i][j] = p
            j += 1
        else:
            i += 1
            matrice_p_interni[i][0] = p
            j = 1

    i = 0
    j = 0

    for p in xyz_esterno:
        if xyz_esterno.index(p) == 0:
            matrice_p_esterni[0][0] = p
            j += 1
        elif j < 20:
            matrice_p_esterni[i][j] = p
            j += 1
            #print(j)
        else:
            i += 1
            matrice_p_esterni[i][0] = p
            j = 1



    # print(matrice_p_interni)

    sottomatrici_interni = np.zeros(((int(len(xyz_interno)/theta_density)-1)*(theta_density-1),4,3))
    i = 0
    for n in range(
            (int(len(xyz_interno)/theta_density)-1)
    ):
        for m in range(theta_density-1):

            sottomatrici_interni[i] = [
                matrice_p_interni[n][m], matrice_p_interni[n][m+1],
                matrice_p_interni[n+1][m], matrice_p_interni[n+1][m+1]
        ]
            i+=1


    sottomatrici_esterni = np.zeros(((int(len(xyz_esterno) / theta_density) - 1) * (theta_density - 1), 4, 3))
    i = 0
    for n in range(
            (int(len(xyz_esterno) / theta_density) - 1)
    ):
        for m in range(theta_density - 1):
            sottomatrici_esterni[i] = [
                matrice_p_esterni[n][m], matrice_p_esterni[n][m + 1],
                matrice_p_esterni[n + 1][m], matrice_p_esterni[n + 1][m + 1]
            ]
            i += 1



    fig = plt.figure()
    ax = Axes3D(fig, auto_add_to_figure=False)
    ax.set_xlim3d(-40, 40)
    ax.set_ylim3d(-40, 40)
    ax.set_zlim3d(0, 40)
    fig.add_axes(ax)

    color = 'w'
    for p in sottomatrici_interni:
        if color == 'w':
            color = 'k'
        else:
            color = 'w'
        verts_temp1 = [[tuple(p[0]), tuple(p[1]), tuple(p[2])]]
        #print(verts_temp1)
        verts_temp2 = [[tuple(p[1]), tuple(p[2]), tuple(p[3])]]
        tri1 = Poly3DCollection(verts_temp1)
        tri2 = Poly3DCollection(verts_temp2)
        tri1.set_color(color)
        tri2.set_color(color)
        ax.add_collection3d(tri1)
        ax.add_collection3d(tri2)
    #verts = [sottomatrici_interni[15]]
    #print(verts)
    #ax.add_collection3d(Poly3DCollection(verts))


    color = 'w'

    for p in sottomatrici_esterni:
        if color == 'w':
            color = 'k'
        else:
            color = 'w'
        verts_temp1 = [[tuple(p[0]), tuple(p[1]), tuple(p[2])]]
        # print(verts_temp1)
        verts_temp2 = [[tuple(p[1]), tuple(p[2]), tuple(p[3])]]
        tri1 = Poly3DCollection(verts_temp1)
        tri2 = Poly3DCollection(verts_temp2)
        tri1.set_color(color)
        tri2.set_color(color)
        ax.add_collection3d(tri1)
        ax.add_collection3d(tri2)
    # verts = [sottomatrici_interni[15]]
    # print(verts)
    # ax.add_collection3d(Poly3DCollection(verts))




    plt.show()


    # step = pi/6
    # for p in xyz:
    #    xyz[xyz.index(p)] = ((p[0]*compute_cos(step), p[1], p[2]*compute_sin(step)))


if __name__ == '__main__':
    costruisci_toro(10, 30, 20, 20)
