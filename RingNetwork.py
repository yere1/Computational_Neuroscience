#   N.B. To run this simualtion properly, your computer has to have the NEURON
# library properly installed. You can find download instructions and more info 
# at https://www.neuron.yale.edu/neuron/download#
#   I have used NEURON 7.5 on OS X Yosemite 10.10.5

import sys
import pylab
sys.path.insert(0, "/Applications/NEURON-7.5/nrn/lib/python") #Be careful to include the proper path
from neuron import h,gui

class SimpleCell(object):
    """ 
        This class implements a SimpleCell object. Each SimpleCell is a 
        neuron soma with an attached narrow and long dendrite.
    """
    def __init__(self):
        """
            This method creates an instance of a SimpleCell object. Each 
            attribute is set by making use of a corresponding method.
        """
        self.x=self.y=self.z=0
        self.set_sections()
        self.set_topology()
        self.set_geometry()
        self.set_biophysics()
        self.syn=self.insert_synapse()
        
    def __repr__(self):
        return "SimpleCell_object"
        
    def set_sections(self):
        """ 
            This method uses the h.Section() function from the NEURON library
            to create two cellular compartments, a soma and a dendrite.
        """
        self.soma=h.Section(name='soma', cell=self)
        self.dend=h.Section(name='dend', cell=self)
    
    def set_topology(self):
        """
            This method uses the connect() method for a NEURON object to 
            establish a biophisical connection between the soma and the
            dendritecreated via set_sections(). 
        """
        self.dend.connect(self.soma(1))
        
    def set_geometry(self):
        """
            This method sets the geometric properties (length and diameter) 
            for the sections defined via set_sections(). It makes use of the
            shape_3D() method, defined below, to define 3D coordinates of 
            cellular compartments.
        """
        self.soma.L=12.6157
        self.soma.diam=12.6157
        self.dend.nseg=101
        self.dend.L=200
        self.dend.diam=1
        self.shape_3D()
    
    def set_biophysics(self):
        """
            This method sets biophysical parameters of cellular compartments, 
            such as axial resistance, channel conductances and membrane 
            capacitances. It also defines the set of channels inserted in the
            membrane: 'pas' for passive leak and 'hh' for Hodgkin-Huxley.
        """
        self.soma.Ra=100
        self.dend.Ra=100
        self.soma.cm=1
        self.dend.cm=1
        
        self.soma.insert('hh')
        self.soma.gnabar_hh=0.12
        self.soma.gkbar_hh=0.036
        self.gl_hh=0.0003
        self.el_hh=-54.3
        
        self.dend.insert('pas')
        self.dend.g_pas=0.001
        self.dend.e_pas=-65

        
    def shape_3D(self):
        """
            This method defines the 3D shape of each SimpleCell object, so 
            that it can later be drawn on a plot. It makes use of the NEURON
            pt3dadd() function defined for cellular compartment objects.
        """
        somalen=self.soma.L
        dendlen=self.dend.L
        
        h.pt3dclear(sec=self.soma)
        h.pt3dadd(0,0,0,self.soma.diam, sec=self.soma)
        h.pt3dadd(somalen,0,0,self.soma.diam, sec=self.soma)
        
        h.pt3dclear(sec=self.dend)
        h.pt3dadd(somalen,0,0,self.dend.diam,sec=self.dend)
        h.pt3dadd(somalen+dendlen,0,0,self.dend.diam, sec=self.dend)
        
        self.all=h.SectionList()
        self.all.wholetree(sec=self.soma)
    
    def insert_synapse(self):
        """
            This method inserts an exponential synapse at the midpoint of the 
            cell dendrite.
        """
        syn=h.ExpSyn(self.dend(0.5))
        syn.tau=2
        return syn
        
    def set_position(self,x,y,z):
        """
            This method can change the absolute 3D location of the network by
            translational motion, leaving all relative positons intact. Besides
            the SimpleCell object, it takes three other argumentseach of which 
            is the change in the respective Cartesian coordinate. 
        """
        for section in self.all:
            for point in range(int(h.n3d(sec=section))):
                h.pt3dchange(point, 
                             x-self.x+h.x3d(point),
                             y-self.y+h.y3d(point),
                             z-self.z+h.z3d(point),
                             h.diam3d(point))
        self.x, self.y, self.z=x, y, z
        
    def z_rotation(self, theta):
        """
            This method rotates the SimpleCell object in 3D space around the 
            z-axis by a user-specified number of degrees.
        """
        theta*=(pylab.pi/180)
        for section in self.all:
            for point in range(int(h.n3d(sec=section))):
                xold=h.x3d(point, sec=section)
                yold=h.y3d(point, sec=section)
                xnew=xold*pylab.cos(theta)-yold*pylab.sin(theta)
                ynew=xold*pylab.sin(theta)+yold*pylab.cos(theta)
                h.pt3dchange(point, xnew, ynew, h.z3d(point, sec=section), 
                             h.diam3d(point, sec=section), sec=section)
                 

class Ring(object):
    """
        This class implements a Ring object. A Ring is a set of SimpleCell
        objects, where the soma midpoint of each SimpleCell projects onto
        the dendrite midpoint of the next SimpleCell. The last SimpleCell
        in the stack projects onto the first.
        Each Ring has three attributes: .network (a list of all SimpleCells 
        that are form the Ring), .extent (the number of cells in the Ring), 
        .connections (the list of NetCon objects (i.e. network connections) 
        between cells in the Ring).
    """
    def __init__(self, extent):
        """
            This method creates an instant of a Ring object. It takes one 
            user-defined argument, 'extent', which is the number of cells 
            in the Ring to be initialized. It sets the values of three 
            attributes, .network, .extent, .connections.
        """
        network=[]
        for index in range(extent): # This 'for loop' creates the cells and establishes their
            cell=SimpleCell()  # 3D positions. The last part isn't conceptually important, because actual 
            theta=index*2*pylab.pi/extent # connection topology is defined separately. However, it provides a nice vizualization.
            cell.z_rotation(theta*180/pylab.pi)
            cell.set_position(50*pylab.cos(theta), 50*pylab.sin(theta), 0)
            network.append(cell)
        
        self.network=network
        self.extent=extent
        self.connections=self.set_connections()
    
    def show(self):
        """
            This method allows to view the Ring network.
        """
        graph=h.PlotShape()
        graph.exec_menu('Show Diam')
    
    def set_connections(self):
        """
            This method establishes network connectivity. Specifically, it 
            associates the voltage at the soma midpoint of a cell as stimulus
            to the synapse of the next cell in the stack. The last cell in the
            sequence projects onto the first cell. It returns 'connections',
            the list of network connections it has created.
        """
        connections=[]
        for index in range(self.extent):
            origin=self.network[index-1].soma
            source=self.network[index-1].soma(0.5)._ref_v
            target=self.network[index].syn
            netcon=h.NetCon(source, target, sec=origin)
            netcon.threshold=-40
            netcon.delay=1
            netcon.weight[0]=0.5
            connections.append(netcon)
        return connections
            
            


def set_currentclamp(cell):
    """
        This function takes a particular cell as an arguemnt and attaches a
        current clamp to the distal end (as defined by 'loc') of its dendrite.
        It returns a NEURON clamp object as the output.
    """
    loc=1
    clamp=h.IClamp(cell.dend(loc))
    clamp.dur=5
    clamp.amp=0.1
    clamp.delay=0
    return clamp

def get_recording(cell):
    """
        This function sets up recording variables and runs the simulation. 
        It does so by first declaring vectors, which will come to contain
        voltage and time data, as HOC vector objects. Then, it uses the 
        corresponding HOC method record() for vector objects to define, what 
        data each vector will come to contain. Then, it runs the simulation and 
        returns those vectors as output.
    """
    Vdendd=h.Vector()
    Vdendm=h.Vector()
    Vdendp=h.Vector()
    t=h.Vector()
    Vdendd.record(cell.dend(1)._ref_v)
    Vdendm.record(cell.dend(0.5)._ref_v)
    Vdendp.record(cell.dend(0)._ref_v)
    t.record(h._ref_t)
    
    tstop=100
    h.tstop=tstop
    h.run()

    return Vdendd, Vdendm, Vdendp, t

def show_output(Vdendd, Vdendm, Vdendp, t, new_fig=False):
    """
        This method allows to visualize results of the simulation. Specifically,
        it plots the voltage at distal, medial and proximal points of the
        dendrite (Vdendd, Vdendm, Vdendp) against time.
    """
    if new_fig:
        pylab.figure(figsize=(8,4))
    pylab.plot(t, Vdendd, color='red')
    pylab.plot(t, Vdendm, color='green')
    pylab.plot(t, Vdendp, color='blue')
    pylab.xlabel('Time, ms', fontsize=18)
    pylab.ylabel('Voltage, mV', fontsize=18)
    pylab.legend(['distal', 'medial', 'proximal'])
            
    
#==============================================================================
#   The following few lines bring together all the classes and functions 
# defined above to run the stimulation. Specifically, the user can choose,
# which cell receives an external current injection and which cell (same or 
# other) will have its membrane voltage recorded. The plot will contain three
# voltage traces, volatge recording from the distal, medial and proximal points
# of the recorded cell's dendrite.
#==============================================================================
ring=Ring(10)
clamped_cell=ring.network[0]
recorded_cell=ring.network[1]
clamp=set_currentclamp(clamped_cell)
Vdendd, Vdendm, Vdendp, t=get_recording(recorded_cell)
show_output(Vdendd, Vdendm, Vdendp, t)