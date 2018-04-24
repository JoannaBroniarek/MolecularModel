import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
        self.t_data = [self.delta_t*i for i in range(self.steps-1)]
        self.e_p = []
        self.e_k = []
        self.e_t = []

    def add_data(self, x, v): #popraw to w zaleznosci jak dziala newton w verlecie
    	try:
            self.x_data.append(x)
            self.v_data.append(v)
    	except:
    	    self.x_data.extend(x)
    	    self.v_data.extend(v)

    def add_energy(self, energy, e_type):
        dic = {"potential":self.e_p, "kinetic":self.e_k, "total":self.e_t}
        dic[e_type].append(energy)

    def run(self): #przyjmuje liczbe klatek
        s = 0
        while s < self.steps-1:
            x, v = self.algorithm(self.x_data, self.v_data, self.delta_t)
            self.add_data(x, v)
            potential = self.potential(self.x_data[-1])
            kinetic = self.kinetic_energy()
            self.add_energy(potential, "potential")
            self.add_energy(kinetic, "kinetic")
            self.add_energy(potential+kinetic, "total")
            s+=1
            #plt.plot(self.e_t, "g-")
            #plt.show()

    def kinetic_energy(self):
        m = 1
        v = self.v_data[-1]
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

    def plot_energy(self):
        plt.plot(self.t_data, self.e_t, "g")
        plt.savefig("plot_energy_newton.png")
        #np.savetxt("x.txt",np.array(self.t_data))
        #np.savetxt("y.txt",np.array(self.e_t))


class Algorithm(object):
    @staticmethod
    def newton(x_data, v_data, delta_t): #energia powinna rosnac- wynika to z bledu
        m = 1
        x1 = x_data[-1] + v_data[-1] * delta_t
        v1 = v_data[-1] + -x_data[-1]/m * delta_t #F=-x oscylator harmoniczny
        return x1, v1

    @staticmethod
    def verlet(x_data, v_data, delta_t):
        m = 1
        if len(x_data)==1:
            newton_x = x_data[:]
            newton_v = v_data
            for s in range(1000):
                x, v = Algorithm.newton(newton_x, newton_v, delta_t/1000)
                newton_x.append(x)
                newton_v.append(v)
            return newton_x[-1], newton_v[-1]
        else:
            x0 = x_data[-2]
            v0 = v_data[-2]
            x1 = x_data[-1]
            v1 = v_data[-1]
            x2 = 2*x1 - x0 + (-x1/m) * ((delta_t)**2)
            v2 = (x1 - x0)/delta_t
            return x2, v2

    @staticmethod
    def velocityverlet(x_data, v_data, delta_t):
        m = 1
        if len(x_data)==1:
            newton_x = x_data
            newton_v = v_data
            for s in range(1000):
                x, v = Algorithm.newton(newton_x, newton_v, delta_t/1000)
                newton_x.append(x)
                newton_v.append(v)
            return newton_x[-1], newton_v[-1]
        else:
            x0 = x_data[-2]
            v0 = v_data[-2]
            x1 = x_data[-1]
            v1 = v_data[-1]
            v2 = v1 + ((-x0/m) + (-x1/m))/2 * delta_t
            x2 = x1 + v2*delta_t + (-x1/2*m) * (delta_t)**2
            return x2, v2


    @staticmethod
    def leapfrog(x_data, v_data, delta_t):
        m = 1
        x1 = x_data[-1]
        v1 = v_data[-1]
        v2 = v1 + (-x1/m)*delta_t
        x2 = x1 + v2*delta_t
        return x2, v2

class Potential(object):
    @staticmethod
    def oscylharm(x):
        k = 1
        potential = k*(x**2/2)
        return potential

class Solver(object):
    def __init__(self, x0, v0, delta_t, steps):
        self.X = [x0]
        self.V = [v0]
        self.d_t = delta_t
        self.steps = steps

    def exactly_oscylator(self, A, fi, t):
        w = 1 #sqrt(k/m)
        #fi = np.arctan(v0/x0)
        #A = x0/ np.sin(t + fi)
        x = A * np.sin(w*t + fi)
        v = A * w * np.cos(w*t + fi)
        return x, v

    def run(self):
        t = 0
        fi = np.arctan(self.V[0]/self.X[0])
        A = self.X[0]/ np.sin(t + fi)
        for i in range(self.steps):
            xx, vv = self.exactly_oscylator(A, fi, t)
            self.X.append(xx)
            self.V.append(vv)
            t = t + self.d_t

    def write(self, filename):
        f = open(str(filename)+'.xyz', 'w')
        nrFrame = 0
        for frame in self.X:    #format xyx, cala animacja -- kilka xyz w jednym pliku
            f.write("1" + "\n") #ruch jednego atomu
            f.write(str(nrFrame)+"\n")
            f.write("Atom1"+ "\t" + str(frame) + "\t" + "0.0" + "\t" + "0.0" + "\n")
            nrFrame += 1
        f.close()

osc = Solver(1.0, 1.0, 0.5, 500)
osc.run()
osc.write("oscylator-solver")


s = Simulation(x0=1.0, v0=1.0, potential=Potential.oscylharm, algorithm=Algorithm.newton, steps=500, delta_t=0.05)
s.run()
s.write("oscylator-newton")
s.write_statistics("statistics-newton")
#s.plot_energy()

### SPRAWDZANIE RZEDU ALGORYTMU-----------------------------
delta_x = [s.x_data[i] - osc.X[i] for i in range(500)]
t_data = [0]
for i in range(499):
    t_data.append(t_data[-1]+0.05)

x_error = [delta_x[i]/0.05**2 for i in range(500)]
plt.plot(x_error, t_data, "g-")

#------------------------------------------------


#s1 = Simulation(x0=1.0, v0=1.0, potential=Potential.oscylharm, algorithm=Algorithm.verlet, steps=500, delta_t=0.01)
#s1.run()
#s1.write("oscylator-verlet")
#s1.write_statistics("statistics-verlet")
#s1.plot

#s2 = Simulation(x0=1.0, v0=1.0, potential=Potential.oscylharm, algorithm=Algorithm.verlet, steps=500, delta_t=0.01)
#s2.run()
#s2.write("oscylator-velocityverlet")
#s2.write_statistics("statistics-velocityverlet")
#s2.plot_energy()

## rzedy algorytmow wedlug rozwiniecia Taylora
## blad = =koncowka taylora
## blad/delta_t do kwadratu albo wyzszych stopni az otzrymamay prosta
