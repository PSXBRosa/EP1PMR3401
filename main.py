import numpy as np
from solver import Mesh
import sys

def eq(key, dr, dtheta, r, k1, k2):
    kernel = np.array([[0,dtheta**2*r*(-dr+2*r)/(4*(dtheta**2*r**2+dr**2)),0],
                        [dr**2/(2*(dtheta**2*r**2+dr**2)),0,dr**2/(2*(dtheta**2*r**2+dr**2))],
                        [0,dtheta**2*r*(dr+2*r)/(4*(dtheta**2*r**2+dr**2)),0]]
                        )

    if key == 21 or key == 20:
        kernel = np.array([[0,dtheta**2*r*(-dr+2*r)/4/(dtheta**2*r**2+dr**2),0],
                            [dr**2/(dtheta**2*r**2+dr**2),0,0],
                            [0,dtheta**2*r*(dr+2*r)/4/(dtheta**2*r**2+dr**2),0]])
    if key == 31 or key == 30:
        kernel = np.array([[0,dtheta**2*k1*r**2*(-dr+2*r)/(dtheta**2*r**2+dr**2)/(-dr*k1+dr*k2+2*r*k1+2*r*k2),0],
                            [dr**2/(2*(dtheta**2*r**2+dr**2)),0,dr**2/(2*(dtheta**2*r**2+dr**2))],
                            [0,dtheta**2*k2*r**2*(dr+2*r)/(dtheta**2*r**2+dr**2)/(-dr*k1+dr*k2+2*r*k1+2*r*k2),0]]
                        )
    if key == 41 or key == 40:
        kernel = np.array([[0,dtheta**2*r*(-dr+2*r)/(4*(dtheta**2*r**2+dr**2)),0],
                        [dr**2*k2/(k1+k2)/(dtheta**2*r**2+dr**2),0,dr**2*k1/(k1+k2)/(dtheta**2*r**2+dr**2)],
                        [0,dtheta**2*r*(dr+2*r)/(4*(dtheta**2*r**2+dr**2)),0]]
                        )
    if key == 51 or key == 50:
        kernel = np.array([[0,dtheta**2*k2*r**2*(-dr+2*r)/(dtheta**2*r**2+dr**2)/(-dr*k1+dr*k2+2*r*k1+2*r*k2),0],
                        [dr**2/(2*(dtheta**2*r**2+dr**2)),0,dr**2/(2*(dtheta**2*r**2+dr**2))],
                        [0,dtheta**2*k1*r**2*(dr+2*r)/(dtheta**2*r**2+dr**2)/(-dr*k1+dr*k2+2*r*k1+2*r*k2),0]]
                        )

    k = k2 if str(key)[-1] == '1' else k1
    scaler = dtheta**2*dr**2*r**2/(2*(dtheta**2*r**2+dr**2))/k
    return kernel, (1,1,1,1), scaler

def eq_r(key, dr, dp, r, k1, k2):
    # 0x > ponto interior
    # 1x > ponto esquerdo
    # 2x > ponto superior
    # 3x > ponto da fronteira esquerda
    # 4x > ponto da fronteira superior
    # 5x > ponto da fronteira direita
    # 6x > ponto direito
    # x0 > ponto do meio A
    # x1 > ponto do meio B

    k = k2 if str(key)[-1] == '1' else k1
    size = (1,1,1,1)
    kernel = -k*np.array([ [0,-1/(2*dr),0],
                           [0, 0,       0],
                           [0, 1/(2*dr),0]])
    if key == 11 or key == 10:
        kernel = -k*np.array([ [0, -3/(2*dr),0],
                               [0,  4/(2*dr),0],
                               [0, -1/(2*dr),0]])
        size = (0,2,1,1)
        
    if key == 31 or key == 30:
        kernel = -k*np.array([ [0,  1/(2*dr),0],
                               [0, -4/(2*dr),0],
                               [0,  3/(2*dr),0]])
        size = (2,0,1,1)
    return kernel, size, None

def eq_p(key, dr, dp, r, k1, k2):
    # 0x > ponto interior
    # 1x > ponto esquerdo
    # 2x > ponto superior
    # 3x > ponto da fronteira esquerda
    # 4x > ponto da fronteira superior
    # 5x > ponto da fronteira direita
    # 6x > ponto direito
    # x0 > ponto do meio A
    # x1 > ponto do meio B

    k = k2 if str(key)[-1] == '1' else k1
    size = (1,1,1,1)
    kernel = -k*np.array([ [0,          0, 0],
                        [-1/(2*dp*r),0,1/(2*dp*r)],
                        [0,          0, 0]])
    if key == 41 or key == 40:
        size = (1,1,2,0)
        kernel = -k*np.array([ [0,          0, 0],
                            [1/(2*dp*r),-4/(2*dp*r),3/(2*dp*r)],
                            [0,          0, 0]])
        
    return kernel, size, None

def eq_k(key, dr, dp, r, k1, k2):
    k = k2 if str(key)[-1] == '1' else k1
    kernel = k*np.array([[0,0,0],[0,1,0],[0,0,0]])
    return kernel, (1,1,1,1), None

if __name__ == '__main__':

    # argumento da linha de comando
    arg1 = sys.argv[-1]

    r1,r2,R1,R2,t,T,ka,kb,v0 = 0.05, 0.08, 0.03, 0.11, np.deg2rad(18), np.deg2rad(40), 5e6, 1e5, 50
    dr,dp = (R2-R1)/80, T/40

    if arg1.lower() == '-y':
        mesh_0 = Mesh(r1,r2,R1,R2,t,T,ka,kb,v0,dr,dp,eq)
        mesh_0.run()
        V = mesh_0.V
        np.savetxt('mesh.txt',V)

    else:
        print('distribuição de tensões carregada pelo arquivo \'mesh.txt\'')
        V = np.loadtxt('mesh.txt')

    mesh_1 = Mesh(r1,r2,R1,R2,t,T,ka,kb,v0,dr,dp,eq_r, initial_grid=V)
    mesh_1.run(w=1, method='jacobi', t_lim=1, immutable=False)
    Jr = mesh_1.V
    np.savetxt('jr.txt',Jr)

    mesh_2 = Mesh(r1,r2,R1,R2,t,T,ka,kb,v0,dr,dp,eq_p, initial_grid=V)
    mesh_2.run(w=1, method='jacobi', t_lim=1, immutable=False)
    Jp = mesh_2.V
    np.savetxt('jp.txt',Jp)

    mesh_3 = Mesh(r1,r2,R1,R2,t,T,ka,kb,v0,dr,dp,eq_k, initial_grid=np.ones_like(V))
    mesh_3.run(w=1, method='jacobi', t_lim=1, immutable=False)
    k = mesh_3.V
    np.savetxt('k.txt',k)

    Z = (-(Jr**2 + Jp**2)/k)
    np.savetxt('qponto.txt',Z)

    mesh_3 = Mesh(r1,r2,R1,R2,t,T,110,500,v0,dr,dp,eq, consts=(30,25))
    mesh_3.run(add=-Z)
    T = mesh_3.V
    np.savetxt('T.txt',T)