import numpy as np
from math import sqrt, pi

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, Aer, execute
from qiskit.providers.aer import QasmSimulator, StatevectorSimulator, UnitarySimulator
from qiskit.providers.aer import Aer
from qiskit.quantum_info.operators import Operator

from qiskit import QuantumCircuit, assemble, Aer
from qiskit.visualization import plot_histogram, plot_bloch_vector
from math import sqrt, pi

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.providers.aer import QasmSimulator
from qiskit.visualization import plot_histogram
import matplotlib.pyplot as plt
from qiskit.quantum_info import Statevector
from qiskit.providers.aer.library import save_statevector
from qiskit.providers.aer import Aer

from qiskit import IBMQ

from qmps.ground_state import Hamiltonian

from scipy.linalg import expm
from scipy.optimize import minimize

from qiskit.providers.aer.noise import NoiseModel
'''
provider = IBMQ.load_account()
backend = provider.get_backend('ibm_osaka')
noise_model = NoiseModel.from_backend(backend)
coupling_map = backend.configuration().coupling_map
basis_gates = noise_model.basis_gates'''

X = np.array([[0, 1], [1, 0]])
Z = np.array([[1, 0], [0, -1]])
I = np.eye(2, 2)
XI = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
IX = np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
ZI = np.array([[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]])
IZ = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]])

def zeros(n):
    a = []
    for i in range(n):
        a.append(0)
    return a

def ei(i,n):
	vi = np.zeros(n)
	vi[i] = 1.0
	return vi[:]

class Loss_cirq:
    def __init__(self, N):
        self.N = N
        self.cb_idx = self.N

        def InitCirc():
            r_qu_n = self.N
            one_lay_gate_num = (2*self.N -1) // 2 ##one layer gate number
            #anci_qu_n = one_lay_gate_num * (cnt + 1)
            anci_qu_n = 0
            total_qubit = r_qu_n + anci_qu_n
            circ = QuantumCircuit(total_qubit, 1)
            return circ

        self.__circ = InitCirc()

    def PrintCirc(self):
        print(self.__circ)

    def U2(self, param, qubits, mode):   #D = 2 tensor
        if mode == 0:
            self.__circ.rz(param[0], qubits[0])
            self.__circ.rx(param[1], qubits[0])
            self.__circ.rz(param[2], qubits[0])
            self.__circ.rz(param[3], qubits[1])
            self.__circ.rx(param[4], qubits[1])
            self.__circ.rz(param[5], qubits[1])
            self.__circ.cx(qubits[0], qubits[1])
            self.__circ.ry(param[6], qubits[0])
            self.__circ.cx(qubits[1], qubits[0])
            self.__circ.ry(param[7], qubits[0])
            self.__circ.rz(param[8], qubits[1])
            self.__circ.cx(qubits[0], qubits[1])
            self.__circ.rz(param[9], qubits[0])
            self.__circ.rx(param[10], qubits[0])
            self.__circ.rz(param[11], qubits[0])
            self.__circ.rz(param[12], qubits[1])
            self.__circ.rx(param[13], qubits[1])
            self.__circ.rz(param[14], qubits[1])
            self.__circ.barrier()
        elif mode == 1:
            self.__circ.rx(param[0], qubits[0])
            self.__circ.rz(param[1], qubits[0])
            self.__circ.rx(param[2], qubits[1])
            self.__circ.rz(param[3], qubits[1])
            self.__circ.cx(qubits[1], qubits[0])
            self.__circ.rx(param[4], qubits[0])
            self.__circ.rz(param[5], qubits[0])    
            self.__circ.rx(param[6], qubits[1])
            self.__circ.rz(param[7], qubits[1])         
            self.__circ.barrier()

    def U2_dagger(self, param, qubits, mode):   #D = 2 tensor
        if mode == 0:      
            self.__circ.rz(-param[14], qubits[1])
            self.__circ.rx(-param[13], qubits[1])
            self.__circ.rz(-param[12], qubits[1])
            self.__circ.rz(-param[11], qubits[0])
            self.__circ.rx(-param[10], qubits[0])
            self.__circ.rz(-param[9], qubits[0])
            self.__circ.cx(qubits[0], qubits[1])
            self.__circ.rz(-param[8], qubits[1])
            self.__circ.ry(-param[7], qubits[0])
            self.__circ.cx(qubits[1], qubits[0])
            self.__circ.ry(-param[6], qubits[0])
            self.__circ.cx(qubits[0], qubits[1])
            self.__circ.rz(-param[5], qubits[1])
            self.__circ.rx(-param[4], qubits[1])
            self.__circ.rz(-param[3], qubits[1])
            self.__circ.rz(-param[2], qubits[0])
            self.__circ.rx(-param[1], qubits[0])
            self.__circ.rz(-param[0], qubits[0])
            self.__circ.barrier()
        elif mode == 1:
            self.__circ.rz(-param[7], qubits[1])
            self.__circ.rx(-param[6], qubits[1])
            self.__circ.rz(-param[5], qubits[0])
            self.__circ.rx(-param[4], qubits[0])
            self.__circ.cx(qubits[1], qubits[0])
            self.__circ.rz(-param[3], qubits[1])
            self.__circ.rx(-param[2], qubits[1])    
            self.__circ.rz(-param[1], qubits[0])
            self.__circ.rx(-param[0], qubits[0])         
            self.__circ.barrier()


    def fmps(self, param):
        for i in range(self.N-1):
            self.U2(param, [i, i+1])

    def fmps_degger(self, param):
        for i in range(self.N-1, 0, -1):
            self.U2_dagger(param, [i-1, i])


    def trotter(self, H, H_1, H_N):
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
            if gate_idxs == gate1list[0]:
                self.__circ.unitary(H_1, 0)
            self.__circ.unitary(H, gate_idxs)

        for gate_idxs in gate2list:
            #print('gate_idxs = ', gate_idxs)
            self.__circ.unitary(H, gate_idxs)
            if gate_idxs == gate2list[len(gate2list)-1]:
                self.__circ.unitary(H_N, self.N-1)
        self.__circ.barrier()

    def add_z(self, idx):
        self.__circ.z(idx)
        self.__circ.barrier()

    def measure_z(self, params):
        self.fmps(params)
        #a.trotter(WW)
        self.add_z(1)
        self.fmps_degger(params)
        p = self.StateVectorRes()
        return p#np.abs(p)

    def getStateVector(self):
        simulator = Aer.get_backend('statevector_simulator')
        #simulator = Aer.get_backend('qasm_simulator')
        result = execute(self.__circ, simulator).result()
        phi = result.get_statevector(self.__circ)
        #phi = phi/np.linalg.norm(phi)
        #print(phi)
        #res = phi[0] #* pow(2, 6) / self.N
        return phi

    def StateVectorRes(self):
        simulator = Aer.get_backend('statevector_simulator')
        #simulator = Aer.get_backend('qasm_simulator')
        result = execute(self.__circ, simulator).result()
        phi = result.get_statevector(self.__circ)
        #phi = phi/np.linalg.norm(phi)
        #print(phi)
        res = phi[0] #* pow(2, 6) / self.N
        return res

    def MeasureRes(self):
        simulator = Aer.get_backend('qasm_simulator')
        shot_num = 8192
        self.__circ.measure_all()
        result = execute(self.__circ, simulator, shots = shot_num).result()

        #if rank == 0:
        target = ''
        for i in range(self.N):
            target = target + '0'
        #print(result.get_counts(self.__circ))
        counts = result.get_counts(self.__circ)['000 0']
        p = (float(counts) / shot_num) # get amplitude
        print('experiment:   ', p)
        return p

    def NoiseMeasureRes(self):
        #provider = IBMQ.load_account()
        #backend = provider.get_backend('ibm_osaka')
        #simulator = Aer.get_backend('qasm_simulator')
        #noise_model = NoiseModel.from_backend(backend)
        #print(noise_model)
        # Get coupling map from backend
        #coupling_map = backend.configuration().coupling_map
        #basis_gates = noise_model.basis_gates
        shot_num = 8192
        self.__circ.measure_all()
        result = execute(self.__circ, Aer.get_backend('qasm_simulator'),
                 coupling_map=coupling_map,
                 basis_gates=basis_gates,
                 noise_model=noise_model, initial_layout = [0, 1, 2], shots = shot_num).result()

        #if rank == 0:
        target = ''
        for i in range(self.N):
            target = target + '0'
        #print(result.get_counts(self.__circ))
        counts = result.get_counts(self.__circ)['000 0']
        p = (float(counts) / shot_num) # get amplitude
        #print('experiment:   ', p)
        return p

    def LocalStateVectorRes(self, idx):
        
        self.__circ.save_probabilities([idx])
        #self.__circ.measure(idx, 0)
        #circuit = transpile(self.__circ)
        #print(circuit)
        
        #backend = Aer.get_backend('statevector_simulator')
        simulator = Aer.get_backend('statevector_simulator')
        snapshot = simulator.run(self.__circ)
        sim_data = snapshot.result().data()
        #print(sim_data['probabilities'])
        return sim_data['probabilities'][0]

    def LocalMeasureRes(self, idx):
        self.__circ.measure(idx, 0)
        self.__circ.save_probabilities([idx])
        #print(self.__circ)
        #circuit = transpile(self.__circ)
        #print(circuit)
        shot_number = 8192
        
        #backend = Aer.get_backend('statevector_simulator')
        simulator = Aer.get_backend('qasm_simulator')
        snapshot = simulator.run(self.__circ, shots = shot_number)
        sim_data = snapshot.result().data()
        #print(sim_data)
        return sim_data['probabilities'][0]

    def OutputCircuit(self):
        return self.__circ

    def MeasureAll(self):
        self.__circ.measure_all()

def measure_overlap(new_params, original_params, mat, mat1, mat2, n = 3, mode = 1):
    l = 15
    if mode == 1:
        l = 8
    L = Loss_cirq(n)
    L.U2(original_params[:l], [0, 1], mode)
    L.U2(original_params[l:], [1, 2], mode)
    L.trotter(mat, mat1, mat2)

    L.U2_dagger(new_params[l:], [1, 2], mode)
    L.U2_dagger(new_params[:l], [0, 1], mode)
    #phi = pow(abs(L.StateVectorRes()), 2)
    #phi = L.NoiseMeasureRes()
    phi = L.MeasureRes()
    #L.PrintCirc()
    #print(-np.abs(phi))
    return -phi

def measure_overlap_circuit(new_params, original_params, mat, mat1, mat2, n = 3, mode = 1):
    l = 15
    if mode == 1:
        l = 8
    L = Loss_cirq(n)
    L.U2(original_params[:l], [0, 1], mode)
    L.U2(original_params[l:], [1, 2], mode)
    L.trotter(mat, mat1, mat2)

    L.U2_dagger(new_params[l:], [1, 2], mode)
    L.U2_dagger(new_params[:l], [0, 1], mode)
    L.MeasureAll()
    return L.OutputCircuit()

def measure_needed_circuit(InputCircuit):
    l = len(InputCircuit)
    simulator = Aer.get_backend('qasm_simulator')
    shot_num = 1024
    result = execute(InputCircuit, simulator, shots = shot_num).result()
    

    value = []
    for i in range(l):
        counts = result.get_counts(InputCircuit[i])['000 0']
        p = (float(counts) / shot_num) # get amplitude
        value.append(-p)
    #print(value)
    return value


def local_measure_overlap(new_params, original_params, mat, mat1, mat2, n = 3, mode = 1):
    l = 15
    if mode == 1:
        l = 8
    local_loss_fn = 0
    for i in range(n):
        L = Loss_cirq(n)
        L.U2(original_params[:l], [0, 1], mode)
        L.U2(original_params[l:], [1, 2], mode)
        L.trotter(mat, mat1, mat2)

        L.U2_dagger(new_params[l:], [1, 2], mode)
        L.U2_dagger(new_params[:l], [0, 1], mode)
        phi = L.LocalMeasureRes(i)
        #L.PrintCirc()
        #print(-np.abs(phi))
        local_loss_fn -= np.abs(phi)/n
    return local_loss_fn
'''
def local_measure_overlap(new_params, original_params, mat, mat1, mat2 n = 3, mode = 1):
    l = 15
    if mode == 1:
        l = 8
    local_loss_fn = 0
    Circuit = []
    for i in range(n):
        L = Loss_cirq(n)
        L.U2(original_params[:l], [0, 1], mode)
        L.U2(original_params[l:], [1, 2], mode)
        L.trotter(mat, mat1, mat2)

        L.U2_dagger(new_params[l:], [1, 2], mode)
        L.U2_dagger(new_params[:l], [0, 1], mode)
        L.MeasureAll()
        #phi = L.LocalMeasureRes(i)
        #L.PrintCirc()
        #print(-np.abs(phi))
        #local_loss_fn -= np.abs(phi)/n
        Circuit.append(L.OutputCircuit())
    value = measure_needed_circuit(Circuit)
    return sum(value)/(l*n)'''

def measure_z_1(param, mode = 1):
    l = 15
    if mode == 1:
        l = 8
    L = Loss_cirq(3)
    L.U2(param[:l], [0, 1], mode)
    L.U2(param[l:], [1, 2], mode)
    phi = np.array(L.getStateVector())[::-1]
    print(phi)
    sigma_z = [[1, 0], [0, -1]]
    sigma_z_matrix = sigma_z
    for i in range(0, 2):
            sigma_z_matrix = np.kron(sigma_z_matrix, np.eye(2, 2))
        #print(sigma_z_matrix)
    phi_sigma = sigma_z_matrix.dot(phi)
    p = phi.conjugate().dot(phi_sigma).real
    return p

def measure_z_middle(param, mode = 1):
    l = 15
    if mode == 1:
        l = 8
    L = Loss_cirq(3)
    L.U2(param[:l], [0, 1], mode)
    L.U2(param[l:], [1, 2], mode)
    phi = np.array(L.getStateVector())[::-1]
    print(phi)
    sigma_z = [[1, 0], [0, -1]]
    sigma_z_matrix = []
    for i in range(0, 3):
        sigma_z_matrix = np.eye(2, 2)
        for j in range(1, 1):
            sigma_z_matrix = np.kron(sigma_z_matrix, np.eye(2, 2))
        sigma_z_matrix = np.kron(sigma_z_matrix, sigma_z)
        for j in range(2, 3):
            sigma_z_matrix = np.kron(sigma_z_matrix, np.eye(2, 2))
        #print(sigma_z_matrix)
    phi_sigma = sigma_z_matrix.dot(phi)
    p = phi.conjugate().dot(phi_sigma).real
    return p

def measure_z_3(param, mode = 1):
    l = 15
    if mode == 1:
        l = 8
    L = Loss_cirq(3)
    L.U2(param[:l], [0, 1], mode)
    L.U2(param[l:], [1, 2], mode)
    phi = np.array(L.getStateVector())[::-1]
    print(phi)
    sigma_z = [[1, 0], [0, -1]]
    sigma_z_matrix = np.eye(2, 2)
    for i in range(0, 1):
            sigma_z_matrix = np.kron(sigma_z_matrix, np.eye(2, 2))
        #print(sigma_z_matrix)
    sigma_z_matrix = np.kron(sigma_z_matrix, sigma_z)
    phi_sigma = sigma_z_matrix.dot(phi)
    p = phi.conjugate().dot(phi_sigma).real
    return p

# This function calculate overlap and gradient of the overlap
def compute_overlap_and_gradient(params, shift, mat, mat1, mat2):
    len_param = len(params)
    E = 0
    E = measure_overlap(params, params + shift, mat, mat1, mat2)
    Circuit = []
    for i in range(len_param):
        Circuit.append(measure_overlap_circuit(params + shift + ei(i,len_param)*np.pi/2.0, params, mat, mat1, mat2))
        Circuit.append(measure_overlap_circuit(params + shift - ei(i,len_param)*np.pi/2.0, params, mat, mat1, mat2))
    
    result = measure_needed_circuit(Circuit)
    #print('E = ', E)
    #print('result = ', result)
    g = []
    for i in range(len_param):
        rplus  = result[2*i]
        rminus = result[1+2*i]
        g.append((rplus-rminus)/2.0)
    #print("g = ", g)
    return E, g

def local_compute_overlap_and_gradient(params, shift, mat, mat1, mat2):
    len_param = len(params)
    E = 0
    result = []
    E = local_measure_overlap(params, params + shift, mat, mat1, mat2)
    for i in range(len_param):
        result.append(local_measure_overlap(params + shift + ei(i,len_param)*np.pi/2.0, params, mat, mat1, mat2))
        result.append(local_measure_overlap(params + shift - ei(i,len_param)*np.pi/2.0, params, mat, mat1, mat2))
    #print('E = ', E)
    #print('result = ', result)
    g = []
    for i in range(len_param):
        rplus  = result[2*i]
        rminus = result[1+2*i]
        g.append((rplus-rminus)/2.0)
    #print("g = ", g)
    return E, g

def compute_overlap_and_gradient_spsa(params,shift,mat, mat1, mat2, count):

    nparameters = len(params)
    # build dictionary of parameters to values
    # {left[0]: parameters[0], .. ., right[0]: parameters[0] + shift[0], ...}

    # Define hyperparameters
    c  = 0.1
    a  = 0.16
    A  = 1
    alpha  = 0.602
    gamma  = 0.101

    a_k = a/np.power(A+count,alpha)
    c_k = c/np.power(count,gamma)

    # Determine the random shift

    delta = np.random.binomial(1,0.5,size=nparameters)
    delta = np.where(delta==0, -1, delta) 
    delta = c_k*delta

    E = 0
    results = []
    E = measure_overlap(params, params + shift, mat, mat1, mat2)
    # Now evaluate the circuits with the parameters assigned
    results.append(measure_overlap(params + shift + delta, params, mat, mat1, mat2))
    results.append(measure_overlap(params + shift - delta, params, mat, mat1, mat2))

    g = np.zeros(nparameters)

    # and the gradient

    rplus  = results[0]
    rminus = results[1]

    for i in range(nparameters):
        # G      = (Ep - Em)/2Δ_i
        # var(G) = var(Ep) * (dG/dEp)**2 + var(Em) * (dG/dEm)**2
        g[i] = a_k*(rplus-rminus)/(2.0*delta[i])

    return E,g 

def local_compute_overlap_and_gradient_spsa(params,shift,mat, mat1, mat2,count):

    nparameters = len(params)
    # build dictionary of parameters to values
    # {left[0]: parameters[0], .. ., right[0]: parameters[0] + shift[0], ...}

    # Define hyperparameters
    c  = 0.1
    a  = 0.16
    A  = 1
    alpha  = 0.602
    gamma  = 0.101

    a_k = a/np.power(A+count,alpha)
    c_k = c/np.power(count,gamma)

    # Determine the random shift

    delta = np.random.binomial(1,0.5,size=nparameters)
    delta = np.where(delta==0, -1, delta) 
    delta = c_k*delta

    E = 0
    results = []
    E = local_measure_overlap(params, params + shift, mat, mat1, mat2)
    # Now evaluate the circuits with the parameters assigned
    results.append(local_measure_overlap(params + shift + delta, params, mat, mat1, mat2))
    results.append(local_measure_overlap(params + shift - delta, params, mat, mat1, mat2))

    g = np.zeros(nparameters)

    # and the gradient

    rplus  = results[0]
    rminus = results[1]

    for i in range(nparameters):
        # G      = (Ep - Em)/2Δ_i
        # var(G) = var(Ep) * (dG/dEp)**2 + var(Em) * (dG/dEm)**2
        g[i] = a_k*(rplus-rminus)/(2.0*delta[i])

    return E,g 

def adam_gradient(count,m,v,g, shift):
    ## This function implements adam optimizer
    beta1 = 0.9
    beta2 = 0.999
    eps   = 1e-8
    alpha = [0.001 for i in range(len(g))]
    if count == 0:
        count = 1

    new_shift = [0 for i in range(len(g))]

    for i in range(len(g)):
        m[i] = beta1 * m[i] + (1 - beta1) * g[i]
        v[i] = beta2 * v[i] + (1 - beta2) * np.power(g[i],2)

        alpha[i] = alpha[i] * np.sqrt(1 - np.power(beta2,count)) / (1 - np.power(beta1,count))

        new_shift[i] = shift[i] - alpha[i]*(m[i]/(np.sqrt(v[i])+eps))

    return new_shift
    

data_path = "/home/tonytsao/tony/thesis/data/"
def fvqte(V, g, h, cost, grad, gd_method, time_steps, max_iter, file_name):
    H0 = Hamiltonian({'ZZ':V, 'X':g, 'Z':h})
    dt = 0.1
    n_steps = time_steps
    WW = expm(-1j*H0.to_matrix()*dt)
    H_1 = np.add(0.5*g*X, 0.5*h*Z)
    W1 = expm(-1j*H_1*dt)
    H_N = np.add(0.5*g*X, 0.5*h*Z)
    WN = expm(-1j*H_N*dt)
    #initial_params = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29]
    initial_params = zeros(16)
    print(initial_params)
    original_params = initial_params
    params = np.array(initial_params)

    shift = np.zeros(16)
    for i in range(16):
        shift[i] = 0.05
    print(shift)

    ths = 0.9999

    #compute_overlap_and_gradient(initial_params, shift, WW)
    evz = 0
    f_evz = open(data_path + file_name, "w")
    evz1 = measure_z_1(params)
    evz = measure_z_middle(params)
    evz3 = measure_z_3(params)
    print('sigma_z = ', evz)
    f_evz.write(str(evz1) + '\n')    
    f_evz.write(str(evz) + '\n')
    f_evz.write(str(evz3) + '\n')
    f_evz.flush()
    for i in range(n_steps):
        print('\n================================== \n')
        print("Time slice:",i+1)
        print("Shift before optimizing this step:",shift)
        print("Initial parameters:", params)
        print('\n================================== \n')
        count = 0
        overlap = [0.01,0]
        g_norm = 1
        #old_grad = np.zeros(len(params)) #momentum
        #g = np.zeros(len(params)) #momentum
        if gd_method == 'Adam':
            m = np.zeros(len(params))
            v = np.zeros(len(params))

        if gd_method == 'momentum':
            old_grad = np.zeros(len(params))
            g = np.zeros(len(params))
        while overlap[0] < ths and count < max_iter:         #ths =0.99, max_iter = 50
            print("Shift optimizing step:",count+1)
            count = count +1 

            if gd_method == 'momentum':
                old_grad = np.asarray(g)
            if cost == 'global':
                if grad == 'param_shift':
                    E,g = compute_overlap_and_gradient(params,shift,WW, W1, WN)
                    
                if grad == 'SPSA':
                    E,g = compute_overlap_and_gradient_spsa(params,shift,WW, W1, WN, count) #spsa

            if cost == 'local':
                if grad == 'param_shift':
                    E,g = local_compute_overlap_and_gradient(params,shift,WW, W1, WN)
                    
                if grad == 'SPSA':
                    E,g = local_compute_overlap_and_gradient_spsa(params,shift,WW, W1, WN, count) #spsa

            print("E = ", E)
            print("g = ", g)
            overlap[0] = np.abs(E)
            
            if gd_method == 'Adam':
                print("\n Adam \n")
                meas_grad = np.asarray(g)
                shift = np.asarray(adam_gradient(count,m,v,meas_grad, shift))

            if gd_method == 'sgd':
                '''
                for i in g:
                    i = 0.5*i'''
                shift = shift - g

            if gd_method == 'momentum':
                
                print("Momentum")
                m_grad = 0.6*np.asarray(g) + 0.4*old_grad
                shift = shift - m_grad

            g_norm = np.linalg.norm(g)
            print("g_norm ", g_norm)

        print('\n---------------------------------- \n')

        print("Shift after optimizing:",shift)
        print("New parameters:"        ,params + shift)
        print("New overlap: "          ,overlap[0])

        params = params + shift
        evz1 = measure_z_1(params)
        evz = measure_z_middle(params)
        evz3 = measure_z_3(params)
        print('sigma_z = ', evz)
        f_evz.write(str(evz1) + '\n')    
        f_evz.write(str(evz) + '\n')
        f_evz.write(str(evz3) + '\n')
        f_evz.flush()

def main(sampler = None):   
    V = -1/2
    g = 0.5
    h = 2
    cost = 'global'
    gradient_method = 'SPSA'
    gd_method = 'sgd'
    time_steps = 40
    max_iter = 800
    txt_name = "global_SPSA_sgd800.txt"
    fvqte(V, g, h, cost, gradient_method, gd_method, time_steps, max_iter,txt_name)


if __name__ == "__main__":
    main()