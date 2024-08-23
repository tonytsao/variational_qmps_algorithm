from qiskit import QuantumCircuit, Aer, execute
from qiskit.quantum_info.operators import Operator
import numpy as np
import os


from qmps.time_evolve_tools import put_env_on_right_site_D4, put_env_on_right_site
from qmps.ground_state import Hamiltonian
from scipy.linalg import expm
import matplotlib.pyplot as plt

from scipy.fftpack import fft, ifft

X = np.array([[0, 1], [1, 0]])
Z = np.array([[1, 0], [0, -1]])
I = np.eye(2, 2)
XI = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
IX = np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
ZI = np.array([[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]])
IZ = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]])
class F_cirq:
    def __init__(self, N, H, cnt):
        self.N = N
        self.H = H
        self.cnt = cnt
        self.cb_idx = self.N

        def InitCirc():
            r_qu_n = self.N
            one_lay_gate_num = (2*self.N -1) // 2 ##one layer gate number
            #anci_qu_n = one_lay_gate_num * (cnt + 1)
            anci_qu_n = 0
            total_qubit = r_qu_n + anci_qu_n
            circ = QuantumCircuit(total_qubit)
            return circ

    
        self.__circ = InitCirc()

    def is_odd(self):
        return True if (self.N)%2 == 1 else False
    
    def Append_U(self, mat, idx):
        self.__circ.unitary(mat, idx)

    def Append1Layer_U(self, mat, mat1, mat2):
        def Gate1List():
            gateN = self.N // 2
            lst = []
            for i in range(0, gateN):
                q1_idx = 2 * i
                lst.append([q1_idx, q1_idx + 1])
            #print(lst)
            return lst
    
        def Gate2List():
            gateN1 = self.N // 2
            gateN = (self.N - 1) // 2
            lst = []
            for i in range(0, gateN):
                q1_idx = 2 * i + 1
                #qa_idx = self.N + gateN1 + i
                lst.append([q1_idx, q1_idx + 1])
            #print(lst)
            return lst
        gate1list = Gate1List()
        gate2list = Gate2List()
        for gate_idxs in gate1list:
            #print('gate_idxs = ', gate_idxs)
            #mat = np.add(mat, XI)
            #print(mat)
            if gate_idxs == gate1list[0]:
                self.__circ.unitary(mat1, 0)
            
            self.__circ.unitary(mat, gate_idxs)

        for gate_idxs in gate2list:
            #print('gate_idxs = ', gate_idxs)
            #mat = np.add(mat, IX)
            self.__circ.unitary(mat, gate_idxs)
            if gate_idxs == gate2list[len(gate2list)-1]:
                self.__circ.unitary(mat2, self.N-1)
            


    def Append1Layer_U_dagger(self, mat):
        def Gate1List():
            gateN = self.N // 2
            lst = []
            for i in range(0, gateN):
                q1_idx = 2 * i
                lst.append([q1_idx, q1_idx + 1])
            #print(lst)
            return lst
    
        def Gate2List():
            gateN1 = self.N // 2
            gateN = (self.N - 1) // 2
            lst = []
            for i in range(0, gateN):
                q1_idx = 2 * i + 1
                #qa_idx = self.N + gateN1 + i
                lst.append([q1_idx, q1_idx + 1])
            #print(lst)
            return lst
        gate1list = Gate1List()
        gate2list = Gate2List()
        for gate_idxs in gate2list:
            #print('gate_idxs = ', gate_idxs)
            self.__circ.unitary(mat, gate_idxs)

        for gate_idxs in gate1list:
            #print('gate_idxs = ', gate_idxs)
            self.__circ.unitary(mat, gate_idxs)


    def Append1Layer(self, mat):
        def Gate1List():
            gateN = self.N // 2
            lst = []
            for i in range(0, gateN):
                q1_idx = 2 * i
                lst.append([q1_idx, q1_idx + 1, self.cb_idx])
                self.cb_idx = self.cb_idx + 1
            #print(lst)
            return lst
    
        def Gate2List():
            gateN1 = self.N // 2
            gateN = (self.N - 1) // 2
            lst = []
            for i in range(0, gateN):
                q1_idx = 2 * i + 1
                qa_idx = self.N + gateN1 + i
                lst.append([q1_idx, q1_idx + 1, self.cb_idx])
                self.cb_idx = self.cb_idx + 1
            #print(lst)
            return lst
        gate1list = Gate1List()
        gate2list = Gate2List()
        for gate_idxs in gate1list:
            #print('gate_idxs = ', gate_idxs)
            self.__circ.unitary(mat, gate_idxs)

        for gate_idxs in gate2list:
            self.__circ.unitary(mat, gate_idxs)

        #self.__circ.reset(2)


    def StateVectorRes(self, cnt):
        simulator = Aer.get_backend('statevector_simulator')
        #simulator = Aer.get_backend('qasm_simulator')
        result = execute(self.__circ, simulator).result()
        phi = result.get_statevector(self.__circ)
        #phi = phi/np.linalg.norm(phi)
        #print(phi)
        res = phi[0] #* pow(2, 6) / self.N
        return res

    def getStateVector(self):
        simulator = Aer.get_backend('statevector_simulator')
        #simulator = Aer.get_backend('qasm_simulator')
        result = execute(self.__circ, simulator).result()
        phi = result.get_statevector(self.__circ)
        #phi = phi/np.linalg.norm(phi)
        #print(phi)
        #res = phi[0] #* pow(2, 6) / self.N
        return phi

    def PrintCirc(self):
        print(self.__circ)


def fqte(N, H, H_1, H_N, cnt, dt = 0.1):
    dt = 0.1
    #print(0.5*H)
    #tmp = put_env_on_left_site_D4(0.5*H)
    #tmp = tmp.reshape(2,2,2,2,2,2)
    #print(tmp[:,:,0,:,:,0].reshape(4,4)-H*0.5)
    #tmp = np.linalg.eig(H)
    #print(tmp[0])
    f_cirq = F_cirq(N, H, cnt)
    
    #print(H)
    e_H = expm(-H*1j*dt)
    e_H_1 = expm(-H_1*1j*dt)
    e_H_N = expm(-H_N*1j*dt)
    #print(e_H)
    e_HU = Operator(e_H)#直接放
    e_HU_1 = Operator(e_H_1)
    e_HU_N = Operator(e_H_N)
    #print("e_HU = ", e_HU)
    #e_H_bar = expm(H*1j*dt)
    #e_HU_bar = Operator(e_H_bar)
    #print("e_HU_bar = ", e_HU_bar)
    #tmp = put_env_on_right_site_D4(0.5*e_H)
    #tmp = tmp.reshape(2,2,2,2,2,2)
    #print(tmp)
    
    sigma_x = [[0, 1], [1, 0]]
    '''
    f_cirq.Append_U(sigma_x, 0)
    f_cirq.Append_U(sigma_x, 1)
    f_cirq.Append_U(sigma_x, 2)
    f_cirq.Append_U(sigma_x, 3)
    f_cirq.Append_U(sigma_x, 4)
    f_cirq.Append_U(sigma_x, 5)
    f_cirq.Append_U(sigma_x, 6)
    f_cirq.Append_U(sigma_x, 7)
    f_cirq.Append_U(sigma_x, 8)
    f_cirq.Append_U(sigma_x, 9)
    f_cirq.Append_U(sigma_x, 10)'''


    for i in range(0, cnt):
        f_cirq.Append1Layer_U(e_HU, e_HU_1, e_HU_N)
    f_cirq.PrintCirc()
    wf = np.array(f_cirq.getStateVector())[::-1]
    return wf
    '''
    sigma_z = [[1, 0], [0, -1]]
    sigma_z_matrix = []
    for i in range(0, N):
        if i == 0:
            sigma_z_matrix = np.eye(2, 2)
        elif i == (N-1)//2:
            sigma_z_matrix = np.kron(sigma_z_matrix, sigma_z)
        else:
            sigma_z_matrix = np.kron(sigma_z_matrix, np.eye(2, 2))
    #print(sigma_z_matrix)
    #print("Magnetization measure = ", e_m)
    #wf_sigma = wf.dot(sigma_z_matrix)
    wf_sigma = sigma_z_matrix.dot(wf)
    #print(wf_sigma)
    p = (wf.conjugate()).dot(wf_sigma).real
    #print("magnetization = ", pow(p, 2))
    return p'''

def cal_sigma_z(wf, N):
    #print(wf)
    sigma_z = [[1, 0], [0, -1]]
    p = []
    for i in range(0, N):
        sigma_z_matrix = []
        if i == 0:
            sigma_z_matrix = sigma_z
            for j in range(i, N-1):
                sigma_z_matrix = np.kron(sigma_z_matrix, np.eye(2, 2))
            #print(sigma_z_matrix)
        else:
            sigma_z_matrix = np.eye(2, 2)
            for j in range(1, i):
                sigma_z_matrix = np.kron(sigma_z_matrix, np.eye(2, 2))
            sigma_z_matrix = np.kron(sigma_z_matrix, sigma_z)
            for j in range(i+1, N):
                sigma_z_matrix = np.kron(sigma_z_matrix, np.eye(2, 2))
            #print(sigma_z_matrix)

        wf_sigma = sigma_z_matrix.dot(wf)
        p.append((wf.conjugate()).dot(wf_sigma).real)
        #print(p)
    return p
    

## main function

def trotter(N, J, hx, hz, cnt, dt = 0.1):
    g = hx
    h = hz
    H0 = Hamiltonian({'ZZ':J, 'X':g, 'Z':h})
    H = H0.to_matrix()
    print(H)
    H_1 = np.add(0.5*g*X, 0.5*h*Z)
    print(H_1)
    H_N = np.add(0.5*g*X, 0.5*h*Z)
    print(H_N)
    #dt = 0.1
    #H_1 = H
    #H_N = H
    es = []
    t = []
    #f_e = open(data_path + "trotter.txt", "w")
    for i in range(0, cnt):
        
        e = cal_sigma_z(fqte(N, H, H_1, H_N, i, dt), N)
        
        #f_e.write(str(e) + '\n')
        #f_e.flush()
        
        es.append(e)
        
        t.append(dt*i)
    #print(es)
    '''
    plt.plot(t, es,"-")
    plt.title('Finite qte, initial state: DDD')
    plt.xlabel('t')
    plt.ylabel('sigma z')
    plt.savefig("old_0.75_0.75.jpg")
    plt.show()'''
    return es

e = trotter(3, -0.5, 0.5, 1.75, 50, 0.1)

data_path = '/home/tonytsao/tony/thesis/data/'
f_e = open(data_path + "a.txt", "w")
for i in range(0, len(e)):
    for j in range(3):
        f_e.write(str(e[i][j]) + '\n')
        f_e.flush()
print(e)