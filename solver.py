import numpy as np
import timeit

np.set_printoptions(threshold=np.inf)
np.set_printoptions(precision=2)
np.set_printoptions(suppress=True)

class Mesh:
    def __init__(self,r1,r2,R1,R2,t,T,ka,kb,v0,dr,dp,eqs,initial_grid=None,consts=None):

        # condições fisicas do problema
        self._R1 = R1
        self._dr = dr
        self._dp = dp
        self._ka = ka
        self._kb = kb

        # margens da estrutura
        self._m = round((R2-R1)/dr + 1)
        self._n = round((T - 0)/dp + 1)
        self._v = round((r1-R1)/dr)
        self._t = round((r2-R1)/dr)
        self._w = round((t - 0)/dp)

        # variaveis
        self.eqs = eqs
        self._v0 = v0
        self._time = 0
        
        # cria a malha
        self._make_grid(initial_grid,consts)

    @property
    def V(self):
        # remove padding
        return self._V[1:-1,1:-1]

    def _classify(self,i,j):

        # 0x > ponto interior
        # 1x > ponto esquerdo
        # 2x > ponto superior
        # 3x > ponto da fronteira esquerda
        # 4x > ponto da fronteira superior
        # 5x > ponto da fronteira direita
        # 6x > ponto direito
        # x0 > ponto do meio A
        # x1 > ponto do meio B

        out = 00
        if (i >= self._v) and (i <= self._t) and (j <= self._w):
            if (i == self._v):
                out = 31
            elif (i == self._t):
                out = 51
            elif (j == self._w):
                out = 41
            else:
                out = 1
        elif (i == 1):
            out = 10
        elif (i == self._m):
            out = 60
        elif (j == self._n):
            out = 20
        return out        

    def _make_grid(self, initial_value = None, consts=None):

        # a malha apresenta um conjunto de pontos em i = 0 e i = n que devem ser imutaveis
        # além disso, existem duas colunas adicionais à esquerda e direita dos pontos que
        # serão iterados. Essas colunas adicionais servem para facilitar o código.

        if initial_value is None:
            s,e = consts if consts is not None else (100,0)
            VL = np.full((1,self._n),s,dtype=float)
            VR = np.full((1,self._n),e,dtype=float)
            VM = np.full((self._m-2,self._n),self._v0,dtype=float)
            V  = np.vstack((VL,VM,VR))
        else:
            V = initial_value
        self._V = np.pad(V, (1,), 'edge')
        self._V_old = np.copy(self._V)

    def _get_kernel(self, i, key):

        # kernels são matrizes 3x3 contendo as constantes que devem ser multiplicados
        # os valores dos pontos ao redor do ponto (i,j) a fim de atender a equação di-
        # ferencial que rege o problema 

        dr,dtheta = self._dr,self._dp
        r = self._R1 + i*dr
        k1,k2 = self._ka, self._kb

        kernel = self.eqs(key, dr, dtheta, r, k1, k2)
        
        return kernel

    def _full_iter(self, w=1.75, method='gauss_seidel', immutable=True, add=None):

        # iteração do algorítmo G-S
        np.copyto(self._V_old,self._V)

        flag = 1 if add is not None else 0
        
        if flag:
            scale_i, scale_j = round(add.shape[0]/self._V.shape[0]), round(add.shape[1]/self._V.shape[1])

        i_s = range(2,self._m) if immutable else range(1,self._m+1)
        for i in i_s:
            for j in range(1,self._n+1):
                
                key = self._classify(i,j)
                kernel,slice_size,scaler = self._get_kernel(i,key)

                if flag:
                    cst = scaler * add[scale_i*(i-1)][scale_j*(j-1)]
                else:
                    cst = 0

                c,b,e,d = slice_size

                if method == 'gauss_seidel':
                    slice  = self._V[i-c:i+b+1,j-e:j+d+1]
                elif method == 'jacobi':
                    slice  = self._V_old[i-c:i+b+1,j-e:j+d+1]
                new_val = np.sum(slice*kernel) + cst
                print(self._V)
                input()
                update = w*new_val + (1-w)*self._V_old[i,j]
                self._V[i,j] = update

                # condição de simetria
                if j == 1:
                    self._V[i,0] = update

    def _get_error(self, threshold, tlim):

        # condição de parada
        error = np.max(np.abs(self._V - self._V_old)/(self._V + 1e-8))
        self._error_in_t = error
        if (self._time == 0 or error > threshold) and self._time < tlim:
            return True
        else:
            return False

    def run(self, w=1.75, method='gauss_seidel',t_lim=np.inf, immutable=True, add=None):
        th = 1e-4
        time = timeit.default_timer()
        while self._get_error(th, t_lim):
            self._full_iter(w=w, method=method, immutable=immutable, add=add)
            print(f'iterações: {self._time} | erro: {self._error_in_t:.5f} | cronômetro: {(timeit.default_timer() - time):.2f}',end='\r')
            self._time += 1
        print(f'\nconcluido em t = {self._time} ({(timeit.default_timer() - time):.2f}s) e erro {self._error_in_t:.4f}')