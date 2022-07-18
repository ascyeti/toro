from sympy import *
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt

def costruisci_toro(r, R, r_density, theta_density):
    x, y = symbols("x y")
    f = sqrt(r**2-(sqrt(x**2+y**2)-R)**2)
    # r è il raggio della circonferenza che ruota attorno all'asse z
    # R è il raggio della circonferenza che passa per tutti i centri delle circonferenze di raggio r


    xyz = np.zeros((r_density, theta_density, 3)) # array dei punti del toro

    for i in range(r_density):
        for j in range(theta_density):

            r_temp = R - r + i * (2 * r) / (r_density - 1)
            theta_temp = 2 * pi * j / (theta_density - 1)

            xp = float(r_temp * cos(theta_temp))
            yp = float(r_temp * sin(theta_temp))
            zp = f.subs([(x, xp), (y, yp)])
            if (not zp.is_real) or zp < 10**(-4): #Messo per bilanciare l'appossimazione
                zp = 0

            xyz[i][j] = [xp, yp, zp] # Poi li sbatto dentro l'array




    numero_sottomatrici = (theta_density-1)*(r_density- 1)
    sottomatrici = np.zeros((numero_sottomatrici, 4, 3))

    i=0
    for n in range(r_density - 1):
        for m in range(theta_density - 1):
            sottomatrici[i] = [
                xyz[n][m], xyz[n][m+1],
                xyz[n+1][m], xyz[n+1][m+1]
            ]
            i += 1

    fig = plt.figure()
    ax = Axes3D(fig, auto_add_to_figure=False)
    ax.set_xlim3d(-40, 40)
    ax.set_ylim3d(-40, 40)
    ax.set_zlim3d(0, 40)
    fig.add_axes(ax)

    color = 'w'
    for p in sottomatrici:
        if color == 'w':
            color = 'k'
        else:
            color = 'w'
        verts_temp1 = [[tuple(p[0]), tuple(p[1]), tuple(p[2])]]
        verts_temp2 = [[tuple(p[1]), tuple(p[2]), tuple(p[3])]]
        tri1 = Poly3DCollection(verts_temp1)
        tri2 = Poly3DCollection(verts_temp2)
        tri1.set_color(color)
        tri2.set_color(color)
        ax.add_collection3d(tri1)
        ax.add_collection3d(tri2)





    plt.show()


if __name__ == '__main__':
    costruisci_toro(10, 30, 20, 20)
