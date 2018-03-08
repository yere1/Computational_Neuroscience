#   This script implements the method for neuron architecture simplification as
# described by Marasco A., Limongiello A., and Migliore M. in their paper 
# "Using Strahler`s analysis to reduce up to 200-fold the run time of realistic
# neuron model" published in Scientific Reports in 2013. Their published code 
# uses the HOC language, which has not been in wide use since1980-ies, so
# implemented the same simplificztion procedure in Python, based on the way it
# is described the Methods section of the paper.
#   I have included a sample file from NeuroMorpho.org which contains the 
# data necessaru to run the code. For more info about the format of the file,
# visit NeuroMorpho.org

from pylab import pi, sqrt

specaxres=250
class AnatomyMeasurements(object):
    """
        This class parses the data in a .txt file from NeuroMorpho.org to 
        get record about the morphology of a nervous cell. The NeuroMorpho
        format iamgines each cell as a collection of cylindrical branches  
        connected end-to-end. Data comes in form of points, belonging to a
        certain branch. The .txt file contains Cartesian coordinates for these points,
        the radius of the branching at that point as well as hierarchical information,
        i.e. what is the "parent" point it is connected to.
        Class attributes are .number (i.e. the id number of the point in the original 
        file), .kind (belongs to a dendrite, axon or soma), .x, .y, .z (Cartesian
        coordinates of the point), .rad (radius of the branching at that point), 
        .upstream (the id number of the upstream point), .downstream (the id numbers
        of its immediate downstream point/s).
    """
    def __init__(self, number, kind, xcoord, ycoord, zcoord, rad, upstream):
        """
            This method accepts relevant values and lays them out in respective
            attributes.
        """
        self.number=int(number)
        self.kind=int(kind)
        self.x=float(xcoord) # Save the x coordinate
        self.y=float(ycoord) # Save the y coordinate
        self.z=float(zcoord) # Save the z coordinate
        self.rad=float(rad) # Save the radius 
        self.upstream=int(upstream) # Save the number of the immediate upstream point
        self.downstream=[] # Compute and save the immediate downstream point/s of the compartment (later)
       
    def __repr__(self):
        return 'Anatomy measurement #'+str(self.number)



class StandardCompartment(object):
    """
        This class, based on the raw anatomy measurements of AnatomyMeasurement,
        creates celluar compartments with certain biophysical properties.
        Class attributes are .number (id of the compartment), .length,
        .radius, .axres (axial resistance), .specaxres (specific axial resistance),
        .parent (id of the parent compartment), .kind (axon, dendrite or soma),
        .nseg (number of segments it should be subdivided to), .strahler 
        (its Strahler number due to branching), .cm (membrane capacitance),
        .gl (leak conductance).
    """
    def __init__(self, number, length, radius, axres, specaxres, parent, kind=-1, nseg=1, strahler=0, cm=0, gl=0):
        self.number=number
        self.length=length
        self.rad=radius
        self.axres=axres
        self.surf=2*pi*length*radius
        self.specaxres=250
        self.parent=parent
        self.kind=kind
        self.strahler=strahler
        self.nseg=nseg
        self.cm=cm
        self.gl=gl
        self.children=[]
        
    def find_strahler(self, compartment_list):
        """
                Find the Strahler number of the compartment given the tree it 
            belongs to as an argument.
                This method makes use of recursion. If the compartment has no children,
            its Strahler number is returned as 1. If it has only one child, the compartment
            gets its child's Strahler number. If it has more than one children, then
            (a) two of the children have the same AND maximal Strahler number, and the compartment 
            gets Strahler number of max+1, (b) the maximal Strahler number among children occurs 
            only once, and the compartment gets Strahler number max.
        """
        if len(self.children)==0:
            return 1
        elif len(self.children)==1:
            child=self.children[0]
            return child.find_strahler(compartment_list)
        else:
            children_strahler=[child.find_strahler(compartment_list) for child in self.children]
            children_strahler.sort()
            if children_strahler[-1]==children_strahler[-2]:
                return children_strahler[-1]+1
            else:
                return children_strahler[-1]
        
    def __repr__(self):
        return 'NeuroMorpho-derived compartment #'+str(self.number)


class RawTree(object):
    """
            This class is effectively a list of the compartments, which together 
            form a single cell. To initilize an instance of the class, the user 
            needs to pass on the name of the .txt file, which contains information
            on the compartments parsed in the NeuroMorpho.org Standard Format. 
            The file needs to be preliminarily cleared off of ANY text that does 
            not relate to a specific compartment, e.g. headers, footnotes, empty 
            lines, etc. Once having extracted the raw data from the file, __init__() 
            uses class-specific methods to compute some relevant parameters, such
            as surface area and the list of children for each compartment.
            Class attributes are .measurements (paresed list of all the info
            in the original raw .txt file), .filename (the name of that file),
            .tree (list of StandardCompartment objects comprising a tree).
            
    """
    def __init__(self, filename):
        """
            This method takes a filename and uses the get_measurements() method
            to parse and save the raw data. Then, it uses add_measurement_children() and
            add_tree_children() to get the children for all compartments. With
            get_soma() it makes a few fixes a few irregularities, which arise 
            because the soma is represented slightly differently from the rest of 
            sections. add_strahler() computes the Strahler number for each 
            unit in the tree. add_cm_gl() defines the membrane capacitance, which
            depends on the Strahler number computed above. add_nseg() defines
            how many subcompartments each compartment should contain.
        """
        self.measurements=self.get_measurements(filename)
        self.add_downstream()
        self.tree,_=self.get_branching(self.measurements[3],2)
        self.get_soma(),
        self.add_tree_children()
        self.add_strahler()
        self.get_item(1).strahler=self.get_item(0).strahler
        self.add_cm_gl()
        self.add_nseg()
        self.filename=filename
        
    def get_measurements(self, filename):
        """
            This method takes up a relevant .txt file and parses its data.
        """
        file=open(filename, 'r')
        raw_text=file.read()
        raw_compartment_list=raw_text.split('\n')
        measurement_list=[]
        for line in raw_compartment_list:
            (_, number, kind, xcoord, ycoord, zcoord, rad, upstream)=line.split(' ')
            measurement=AnatomyMeasurements(number, kind, xcoord, ycoord, zcoord, rad, upstream)
            measurement_list.append(measurement)
       
        return measurement_list
   
    def add_downstream(self):
        """
            This method finds the immediate downstream point/s to a given AnatomyMeasurement point.
        """
        for measurement in self.measurements:
            if measurement.kind==1: #This makes sure the soma is treated separately
                continue
            else:
                measurement_upstream=self.measurements[measurement.upstream-1]
                measurement_upstream.downstream.append(measurement)
                
    def get_branching(self, point, tracker=0, parent=0, tree=[]):
        """
            This method finds all the points that belong to the same branch and 
            uses info about them, saved in the AnatomyMeasurement object, to 
            construct a cellular compartment. The compartment's length is defined
            as the length of the polyline traced by its points. The radius of the
            compartment is defined as a weighted average of radii at each measured 
            point.
        """
        l_tot=0
        r_num_tot=0
        while True:
            point_upstream=self.measurements[point.upstream-1]
            l_loc=((point.x-point_upstream.x)**2+(point.y-point_upstream.y)**2+(point.z-point_upstream.z)**2)**0.5
            l_tot+=l_loc
            r_num_tot+=0.5*(point_upstream.rad+point.rad)*l_loc
            if len(point.downstream)==1:
                point=point.downstream[0]
            elif len(point.downstream)==0:
                r_tot=r_num_tot/l_tot
                axres=0.01*(specaxres*l_tot)/(pi*r_tot**2)
                segment=StandardCompartment(tracker, l_tot, r_tot, axres, specaxres, parent, point.kind)
                tree.append(segment)
                break
            elif len(point.downstream)>1:
                r_tot=r_num_tot/l_tot
                axres=0.01*(specaxres*l_tot)/(pi*r_tot**2)
                segment=StandardCompartment(tracker, l_tot, r_tot, axres, specaxres, parent, point.kind)
                tree.append(segment)
                parent=tracker
                for child in point.downstream:
                    tracker+=1
                    (_,tracker)=self.get_branching(child,tracker,parent,tree)
                break

        return tree, tracker
    
            
    def add_tree_children(self):
        """
            After having defined distinct branches, this method finds the
            children branches of a given compartment.
        """
        for compartment in self.tree:
             parent=self.get_item(compartment.parent)
             parent.children.append(compartment)
        self.tree[0].children.remove(self.tree[0])
                
    def add_strahler(self):
        """
            This method iterates over all compartments in the tree to find its
            branches' Strahler numbers.
        """
        compartment_list=self.tree.copy()
        for comp in compartment_list:
            comp.strahler=comp.find_strahler(compartment_list) 
    
    def add_cm_gl(self):
        """
            Once the Strahler number for each compartment is known, this method
            computes the membrane capacitances.
        """
        for comp in self.tree[0:2]:
            comp.cm=0.8
            comp.gl=0.00006
            
        for comp in self.tree[2:]:
            comp.gl=0.00021
            if comp.strahler>3:
                comp.cm=0.8
            else:
                comp.cm=1.5
                
    def add_nseg(self):
        """
            This method determines how many subsections each compartment should have.
        """
        for comp in self.tree:
            lambda_d=0.1
            lambda_f=1e5*sqrt(2)*sqrt(2*comp.rad/(4*pi*100*comp.specaxres*comp.cm))
            nseg=1+2*int(0.5*(comp.length/(lambda_f*lambda_d)+0.9))
            comp.nseg=nseg
    
    def get_soma(self):
        """
            The representation of the soma in NeuroMorpho is slightly different
            from the cell processes. This method helps to resolve those
            inconsistencies by treating the separately.
        """
        comp_parent=self.measurements[0]
        for index in range(0,2):
            comp=self.measurements[index+1]
            l_loc=((comp.x-comp_parent.x)**2+(comp.y-comp_parent.y)**2+(comp.z-comp_parent.z)**2)**0.5
            r_loc=comp.rad
            axres=0.01*(specaxres*l_loc)/(pi*r_loc**2)
            parent=0
            soma=StandardCompartment(index, l_loc, r_loc, axres, specaxres, parent, comp.kind)
            self.tree.insert(index,soma)

    def get_item(self, number):
        """
            A simple getter function for obtaining the StandardCompartment object
            with a given id number. 
        """
        for item in self.tree:
            if item.number==number:
                return item
      


class EquivalentTree():
    """
        This class helps to carry out the simplification of the original cell.
        Based on Strahler numbers of compartments, it groups them into three different
        categories and then either leaves them intact or carries out simplification 
        procedures described in the original paper. It has only one attribute,
        '.tree', which is a list of all StandardCompartment objects comprising
        the original cell.
    """
    def __init__(self, RawTree_object, thresh):
        """
            To create an EquivelentTree with simplified architecture, this method
            needs the original tree in form of the RawTree object, and the 
            threshold Strahler number for simplification.
        """
        raw_tree=RawTree_object.tree
        trunk=[]
        spiny_strahler=range(1,min(4,thresh))
        smooth_strahler=range(4,min(raw_tree[0].strahler, thresh))
        
        tracker=0
        tracking={}
        for comp in raw_tree:
            if comp.strahler>=thresh:
                tracking.update({comp.number:tracker})
                tracker+=1
                
        
        soma_first=raw_tree[0]
        soma_second=raw_tree[1]
        soma=StandardCompartment(1, soma_first.length+soma_second.length, soma_first.rad, soma_first.axres+soma_second.axres, soma_first.specaxres, 0)
        trunk.append(soma)
        
        for comp in raw_tree[2:]:
            if comp.strahler>=thresh:
                
                parent=tracking[comp.parent]
                eq_comp=StandardCompartment(tracking[comp.number], comp.length, comp.rad, comp.axres, comp.specaxres, parent)
                trunk.append(eq_comp)
                parent=tracking[comp.number]
                
                (spiny_leq, spiny_radeq, spiny_axreseq, spiny_specaxreseq)=find_equivalent(comp, spiny_strahler)
                if spiny_leq!=0:
                    eq_spiny=StandardCompartment(tracker, spiny_leq, spiny_radeq, spiny_axreseq, spiny_specaxreseq, parent)
                    trunk.append(eq_spiny)
                    tracker+=1
                    
                
                (smooth_leq, smooth_radeq, smooth_axreseq, smooth_specaxreseq)=find_equivalent(comp, smooth_strahler)
                if smooth_leq!=0:
                    eq_smooth=StandardCompartment(tracker, smooth_leq, smooth_radeq, smooth_axreseq, smooth_specaxreseq, parent)
                    trunk.append(eq_smooth)
                    tracker+=1
                    
                    buds=find_spiny_on_smooth(comp, spiny_strahler, smooth_strahler)
                    eq_spiny_on_smooth=find_equivalent_spiny_bundle(buds, spiny_strahler, tracker)
                    trunk.append(eq_spiny_on_smooth)
                    tracker+=1    
                        
        self.tree=trunk
   
def find_spiny_on_smooth(comp, spiny_strahler, smooth_strahler):
    """
            Given a graph tree, we can imagine two categories of nodes (i.e. compartments), based on their Strahler numbers. Nodes with 
        larger Strahler numbers form "branches", while the smaller ones form "leaves". This function takes up a node of the "branch"
        type, and returns the set of all "leaves" down the road from from the argument PLUS the argument itself, which are root to a 
        "leaf" cluster, while ignoring all the branches down the road that do not give rise to "leaves" OR are not a root of a leaf cluster. 
        If given the root of the "branch", this function can effectively locate all the "leaf" clusters that come out of it, by giving the "branch" 
        node, on which each "leaf" cluster "sits on". This information can be used later to collapse each "leaf" cluster into one equivalent 
        compartment. For our purposes, "branch" is a compartment corresponding to smooth dendrites with Strahler number below the collapsing 
        threshold in range(4,thresh), while "leaf" is a compartment corresponding to spiny dendrites with Strahler number in range(1,4).
    """
    leaf_root=[]
    for child in comp.children:
        if child.strahler in spiny_strahler:
            leaf_root+=[child]
        elif child.strahler in smooth_strahler:
            leaf_root+=find_spiny_on_smooth(child, spiny_strahler, smooth_strahler)
            
    return set(leaf_root)

    
def find_equivalent(comp, relevant):
    """
            This function takes a compartment and replaces all of its relevant progeny (NOT including the compartment itself)
        with a single equiavent compartment. The relevance of progeny depends on its Strahler number: a compartmmet is deemed
        relevant if its Strahler number lies within a certain range. If the argument has no relevant children, then the 
        geometric parameters for the equivalent progeny are returned as zero. If the compartment has only one relevant child,
        then the sequential merging rules apply to this child PLUS, recursively, all of the child's own relevant progeny (i.e. the 
        grandchildren of the original argument). If there are more than one relevant children, then branchpoint merging rules 
        apply to these children PLUS, recursively, the progeny thereof (i.e. the grandchildren of the original argument).
            Inputs: comp /an object with two attributes among others: .children (list of its children) and .strahler (Strahler number)
                    relevant /an iterable containing relevant Strahler numbers
            Outputs: leq /length of the newly formed equivalent compartment
                     radeq /radius of the newly formed equivalent compartment
                     axreseq /axial resistance of the newly formed equivalent compartment
                     specaxreseq /specific axial resistance of the newly formed equivalent compartment
    """
    relevant_children=[child for child in comp.children if child.strahler in relevant]
    if len(relevant_children)==0:
        return (0,0,0,comp.specaxres)
    elif len(relevant_children)==1:
        child=relevant_children[0]
        (leq_recur, radeq_recur, axreseq_recur, specaxreseq_recur)=find_equivalent(child, relevant)
        leq=child.length+leq_recur
        axreseq=child.axres+axreseq_recur
        specaxreseq=(child.specaxres+specaxreseq_recur)/2
        radeq=0.1*sqrt(specaxreseq*leq/(pi*axreseq))
    else:
        leq_num=0
        leq_denom=0
        axres_sum=0
        axres_prod=1
        radeq_sq=0
        specaxreseq=0
        for child in relevant_children:
            (leq_recur, radeq_recur, axreseq_recur, specaxreseq_recur)=find_equivalent(child, relevant)
            surf=pi*(leq_recur+child.length)*2*(0.1*sqrt(0.5*(child.specaxres+specaxreseq_recur)*(child.length+leq_recur)/(pi*(axreseq_recur+child.axres))))
            leq_num+=surf*(leq_recur+child.length)
            leq_denom+=surf
            axres_sum+=child.axres+axreseq_recur
            axres_prod*=child.axres+axreseq_recur
            radeq_sq+=(0.1*sqrt(0.5*(child.specaxres+specaxreseq_recur)*(child.length+leq_recur)/(pi*(axreseq_recur+child.axres))))**2
        leq=leq_num/leq_denom
        axreseq=axres_sum/len(relevant_children)
        radeq=radeq_sq**0.5
        specaxreseq=100*pi*(radeq**2)*axreseq/leq
        
    return leq, radeq, axreseq, specaxreseq


def find_equivalent_spiny_bundle(exp, spiny_strahler, tracker):
    leq_tot=0
    surf_tot=0
    rad_tot=0
    axreseq_tot=0
    for leaf in exp:
        (spiny_on_smooth_leq, spiny_on_smooth_radeq, spiny_on_smooth_axreseq, spiny_on_smooth_specaxreseq)=find_equivalent(leaf, spiny_strahler)
        leq=leaf.length+spiny_on_smooth_leq
        axreseq=leaf.axres+spiny_on_smooth_axreseq
        specaxreseq=(leaf.specaxres+spiny_on_smooth_specaxreseq)/2
        radeq=0.1*sqrt(specaxreseq*leq/(pi*axreseq))
        
        surf=pi*2*radeq*leq
        leq_tot+=surf*leq
        surf_tot+=surf
        rad_tot+=radeq**2
        axreseq_tot+=axreseq
        
    leq_tot/=surf_tot
    rad_tot**=0.5
    axreseq_tot/=len(exp)
    specaxreseq_tot=100*pi*(rad_tot**2)*axreseq_tot/leq_tot
    equivalent_spiny=StandardCompartment(tracker, leq_tot, rad_tot, axreseq_tot, specaxreseq_tot, tracker-1)
    
    return equivalent_spiny


#==============================================================================
#                   Main body of the code
#==============================================================================
OriginalCell=RawTree("UseForStrahlerSimplification.txt")
ReducedCell=EquivalentTree(OriginalCell, 4)
