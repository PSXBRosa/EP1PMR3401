from matplotlib import projections
import numpy as np
import matplotlib.pyplot as plt

def main():

    Jr = np.loadtxt('jr.txt')
    Jp = np.loadtxt('jp.txt')
    Z  = np.loadtxt('qponto.txt')
    

    Jr = np.hstack((np.flip(Jr,axis=1),Jr))
    Jp = np.hstack((np.flip(Jp,axis=1),Jp))
    Z  = np.hstack((np.flip(Z,axis=1),Z))

    m = len(Jr)
    n = len(Jr[0])

    r,p = np.linspace(0.03,0.11,m), np.linspace(-40*np.pi/180,40*np.pi/180,n)
    R,P = np.meshgrid(r,p)
    x = R*np.cos(P)
    y = R*np.sin(P)
    Jx = Jr.T*np.cos(P) - Jp.T*np.sin(P)
    Jy = Jr.T*np.sin(P) + Jp.T*np.cos(P)

    fig = plt.figure(figsize=(8,6))
    ax0 = fig.add_subplot(1,2,1)
    ax0.quiver(x,y,Jx,Jy)
    ax = fig.add_subplot(1,2,2,projection='3d')
    ax.contour3D(x,y,Z.T,200)
    plt.savefig('J e Qponto')
    # ax.plot_surface(x,y,-Z.T)

    plt.show()

if __name__ == '__main__': main()