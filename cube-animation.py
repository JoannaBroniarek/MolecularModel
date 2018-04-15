import numpy as np

class Generator(object):
    @staticmethod
    def scale(factory, index, cube):
        outFrame = cube.animation[index] * factory
        return outFrame

    @staticmethod
    def rotate(angle, index, cube):
        m = np.array([[np.cos(angle),-np.sin(angle),0], [np.sin(angle),np.cos(angle), 0], [0, 0, 1]])
        outFrame = np.dot(cube.animation[index], m)
        return outFrame


class Cube(object):
    def __init__(self, intEdge, lenght):
        self.intEdge = intEdge #liczba atomow na krawedzi
        self.lenght = lenght #odleglosc miedzy atomami
        self.animation = []
        #self.animation = ['<tutaj umieszczam klatki>'] #format klatki: [<x>,<y>,<z>]

    def create(self):
    	a = self.intEdge * self.lenght
    	C = np.linspace(0, a, self.intEdge, endpoint = True)
    	coords = np.array([np.array([x,y,z]) for x in C for y in C for z in C])
    	self.animation.append(coords)	#pierwsza klatka w animacji

    def writeFrame(self,filename):
        frame = self.animation[-1]
        f = open(str(filename)+'.xyz', 'w')
        f.write(str(self.intEdge**3)+"\n\n")
    	count = 0
    	for coord in frame:
    		f.write("Atom"+str(count)+ "\t" + str(coord[0]) + "\t" + str(coord[1]) + "\t" + str(coord[2]) + "\n")
    		count+=1

    def writeAnimation(self,filename):
        f = open(str(filename)+'.xyz', 'w')
        nrFrame = 0
        for frame in self.animation:    #format xyx, cala animacja -- kilka xyz w jednym pliku
            f.write(str(self.intEdge**3) + "\n")
            f.write(str(nrFrame)+"\n")
            nrAtom = 0
            for coord in frame:
                f.write("Atom"+str(nrAtom)+ "\t" + str(coord[0]) + "\t" + str(coord[1]) + "\t" + str(coord[2]) + "\n")
                nrAtom += 1
            nrFrame += 1

    def makeAnimation(self, nMax, generator, *argsGen):
        for n in range(nMax):
            self.animation.append(generator(*argsGen, cube=self))

c = Cube(6, 1)
c.create()
c.makeAnimation(100,Generator.rotate,np.pi/36,-1)
#c.makeAnimation(100,Generator.scale,1+1./100,-1)
c.writeAnimation("testowaanimacja")
