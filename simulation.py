import numpy as np
#from scipy.misc import derivative

class Simulation(object): #dla jednego atomu
    def __init__(self, x0, v0, potential, algorithm, steps, delta_t): #potencjal i algorytm to funkcje
        self.x0 = x0
        self.v0 = v0
        self.potential = potential
        self.algorithm = algorithm
        self.steps = steps
        self.delta_t = delta_t
        self.x_data = [self.x0]
        self.v_data = [self.v0]
        self.e_p = []
        self.e_k = []
        self.e_t = []

    def add_data(self, x, v):
        self.x_data.append(x)
        self.v_data.append(v)

    def add_energy(self, energy, e_type):
        dic = {"potential":self.e_p, "kinetic":self.e_k, "total":self.e_t}
        dic[e_type].append(energy)

    def run(self): #przyjmuje liczbe klatek
        s = 0
        while s < self.steps:
            x, v = self.algorithm(self.x_data, self.v_data, self.delta_t)
            self.add_data(x, v)
            potential = self.potential(x)
            kinetic = Simulation.kinetic_energy(v)
            self.add_energy(potential, "potential")
            self.add_energy(kinetic, "kinetic")
            self.add_energy(potential+kinetic, "total")
            s+=1

    @staticmethod
    def kinetic_energy(v):
        m = 1
        energy = (m * (v**2))/2
        return energy

    def write_statistics(self, filename):#liczy energie calk dla kolejnych krokow czasowych
        f = open(str(filename)+'.csv', 'w')
        f.write("time\tkinetic\tpotential\ttotal\n")
        time = 0
        for i in range(len(self.x_data)-1):
            f.write(str(time)+"\t"+str(self.e_k[i])+"\t"+str(self.e_p[i])+"\t"+str(self.e_t[i])+"\n")
            time+=self.delta_t
        f.close()

    def write(self, filename):
        f = open(str(filename)+'.xyz', 'w')
        nrFrame = 0
        for frame in self.x_data:    #format xyx, cala animacja -- kilka xyz w jednym pliku
            f.write("1" + "\n") #ruch jednego atomu
            f.write(str(nrFrame)+"\n")
            f.write("Atom1"+ "\t" + str(frame) + "\t" + "0.0" + "\t" + "0.0" + "\n")
            nrFrame += 1
        f.close()

class Algorithm(object):
    @staticmethod
    def newton(x_data, v_data, delta_t):
        m = 1
        x = x_data[-1] + v_data[-1] * delta_t
        v = v_data[-1] + -x/m * delta_t #F=-x oscylator harmoniczny
        #d_x = np.abs(x_data[-1] - x)
        #v = v_data[-1] + derivative(potential, x_data[-1]) * delta_t
        #print x, v
        return x, v

    @classmethod
    def verlet(cls, delta_t):
        count = 0
        m = 1 #masa umownie rowna 1
        x1 = cls.x_data[-1]
        #pierwszy krok dzielimy na 1000 czesci i obliczamy newtonem
        if count==0:
            new_t = delta_t/1000
            x1 = x0 + v0 * new_t
            v1 = v0 + a0 * new_t
            cls.update(x1, v1)
        x = 2*cls.x_data[-1] - x0 + a0 * (delta_t)**2
        v = (x2 - x1)/delta_t
        return x, v

    @staticmethod
    def velocityverlet(x0, v0, delta_t):
        pass

    @staticmethod
    def leapfrog(x0, v0, delta_t):
        pass

class Potential(object):
    @staticmethod
    def oscylharm(x):
        k = 1
        potential = k*(x**2/2)
        return potential


s = Simulation(x0=1.0, v0=1.0, potential=Potential.oscylharm, algorithm=Algorithm.newton, steps=100, delta_t=0.1)
s.run()
s.write("oscylatorharmoniczny")
s.write_statistics("statistics")
