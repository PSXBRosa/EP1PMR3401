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

    return kernel, (1,1,1,1), ""

def eq_temp(key, dr, dtheta, r, k1, k2):
    k = k2 if str(key)[-1] == '1' else k1
    scaler = dtheta**2*dr**2*r**2/(2*(dtheta**2*r**2+dr**2))/k

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

    if key == 10 or key == 11:
        kernel = np.array([[0,0,0],[0,1,0],[0,0,0]])

    if key == 60 or key == 61:
        h = 50
        denominator = dtheta**2*dr**2*h*r + 2*dtheta**2*dr*h*r**2 + 2*dtheta**2*k*r**2 + 2*dr**2*k
        kernel = np.array([[0,2*dtheta**2*k*r**2/denominator,0],
                          [dr**2*k/denominator,0,dr**2*k/denominator],
                          [0,dtheta**2*dr*h*r*(dr+2*r)/denominator,0]])
        scaler =  dtheta**2*dr**2*r**2/denominator

    if key == 600 or key == 610:
        h = 50
        denominator = dtheta**2*dr**2*h*r + 2*dtheta**2*dr*h*r**2 + 2*dtheta**2*k*r**2 + 2*dr**2*k
        kernel = np.array([[0,2*dtheta**2*k*r**2/denominator,0],
                          [2*dr**2*k/denominator,0,0],
                          [0,dtheta**2*dr*h*r*(dr+2*r)/denominator,0]])
        scaler =  dtheta**2*dr**2*r**2/denominator

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
        
    if key == 61 or key == 60:
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
    if key == 21 or key == 20:
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

    r1,r2,R1,R2,t,T,ka,kb,v0 = 0.05, 0.08, 0.03, 0.11, np.deg2rad(18), np.deg2rad(40), 5e-6, 1e-5, 50
    dr,dp = (R2-R1)/80, T/48

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

    mesh_3 = Mesh(r1,r2,R1,R2,t,T,110,500,25+273,2*dr,2*dp,eq_temp, consts=(30+273,25+273))
    mesh_3.run(add = -Z, immutable=False)
    Temp = mesh_3.V
    np.savetxt('T.txt',Temp)

    mesh_4 = Mesh(r1,r2,R1,R2,t,T,110,500,0,2*dr,2*dp,eq_r, initial_grid=Temp)
    mesh_4.run(w=1, method='jacobi', t_lim=1, immutable=False)
    Qr = mesh_4.V
    np.savetxt('Qr.txt',Qr)

    mesh_5 = Mesh(r1,r2,R1,R2,t,T,110,500,0,2*dr,2*dp,eq_p, initial_grid=Temp)
    mesh_5.run(w=1, method='jacobi', t_lim=1, immutable=False)
    Qp = mesh_5.V
    np.savetxt('Qp.txt',Qp)

    # J na direção radial em r = Rx: Jr_Rx
    Jr_R1 = Jr[ 0,:]
    Jr_R2 = Jr[-1,:]
    Qr_R2 = Qr[-1,:]

    # integral pela regra do trapézio composta com h = R1*dp
    I_R1 = (np.sum(Jr_R1) - Jr_R1[0]/2 - Jr_R1[-1]/2)*R1*dp

    # erro estimado da integração
    erro_I_R1 = -1/12*(R1*dp)**3*(len(Jr_R1)-1)*np.min([(Jr_R1[2] - 2*Jr_R1[1] + Jr_R1[0])/(R1*dp)] + [(Jr_R1[i+1] - 2*Jr_R1[i] + Jr_R1[i-1])/(R1*dp) for i in range(1,len(Jr_R1) - 1)] + [(Jr_R1[-1] - 2*Jr_R1[-2] + Jr_R1[-3])/(R1*dp)])

    # integral pela regra do trapézio composta com h = R2*dp
    I_R2 = (np.sum(Jr_R2) - Jr_R2[0]/2 - Jr_R2[-1]/2)*R2*dp

    # erro estimado da integração
    erro_I_R2 = -1/12*(R2*dp)**3*(len(Jr_R2)-1)*np.min([(Jr_R2[2] - 2*Jr_R2[1] + Jr_R2[0])/(R2*dp)] + [(Jr_R2[i+1] - 2*Jr_R2[i] + Jr_R2[i-1])/(R2*dp) for i in range(1,len(Jr_R2) - 1)] + [(Jr_R2[-1] - 2*Jr_R2[-2] + Jr_R2[-3])/(R2*dp)])

    print(f'Corrente I calculada em r = R1: {I_R1:.2e} com erro estimado de {erro_I_R1:.3e}')
    print(f'Corrente I calculada em r = R2; {I_R2:.2e} com erro estimado de {erro_I_R2:.3e}')

    print(f'Resistência média avaliada com I calculado em r = R1: {100/I_R1} e com r = R2: {100/I_R2}')

    # integral pela regra do trapézio composta com h = R2*dp
    q_conv = (np.sum(Qr_R2) - Qr_R2[0]/2 - Qr_R2[-1]/2)*R2*dp

    # erro estimado da integração
    erro_q = -1/12*(R2*dp)**3*(len(Qr_R2)-1)*np.min([(Qr_R2[2] - 2*Qr_R2[1] + Qr_R2[0])/(R2*dp)] + [(Qr_R2[i+1] - 2*Qr_R2[i] + Qr_R2[i-1])/(R2*dp) for i in range(1,len(Qr_R2) - 1)] + [(Qr_R2[-1] - 2*Qr_R2[-2] + Qr_R2[-3])/(R2*dp)])

    print(f'quantidade de calor pela parade de convecção: {q_conv:.2e} com erro estimado de {erro_q:.3e}')




    