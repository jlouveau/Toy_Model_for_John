# This is a Python port of Joy's affinity maturation code
# available at: https://github.com/jlouveau/Toy_Model_for_John

import sys
import numpy as np                          # numerical tools
from copy import deepcopy                   # deepcopy copies a data structure without any implicit references
from scipy.stats import genextreme          # generalized extreme value distribution
from timeit import default_timer as timer   # timer for performance


###### Global parameters ######
    
p_mut        = 0.2                              # probability of mutation per division round
p_CDR        = 1.0                              # probability of mutation in the CDR region
p_CDR_lethal = 0.3                              # probability that a CDR mutation is lethal
p_CDR_silent = 0.5                              # probability that a CDR mutation is silent
p_CDR_affect = 1. - p_CDR_lethal - p_CDR_silent # probability that a CDR mutation affects affinity
p_FR_lethal  = 0.8                              # probability that a framework (FR) mutation is lethal
p_FR_silent  = 0.                               # probability that a FR mutation is silent
p_FR_affect  = 1. - p_FR_lethal - p_FR_silent   # probability that a FR mutation affects affinity

nb_Ag        = 1               # number of antigens
conc         = 1.20            # antigen concentration
energy_scale = 0.30            # inverse temperature
help_cutoff  = 0.70            # only B cells in the top (help_cutoff) fraction of binders receive T cell help
p_recycle    = 0.70            # probability that a B cell is recycled
p_exit       = 1. - p_recycle  # probability that a B cell exits the GC


# CHECK FOR CONGRUENCE WITH MATLAB
c      = 0.7                                    # Generalized extreme value (GEV) distribution shape parameter
scale  = 1.2                                    # GEV scale parameter
loc    = -1.5                                   # GEV location parameter
gevrnd = genextreme(c=c, scale=scale, loc=loc)  # GEV random number generator

mu     = 1.9    # lognormal mean
sigma  = 0.5    # lognormal standard deviation
o      = 3.0    # lognormal offset


###### B cell clone class ######

class BCell:

    def __init__(self, nb = 512, **kwargs):
        """ Initialize clone-specific variables. """
        self.nb = nb                                        # default starting population size = 512 (9 rounds of division)
        
        if ('Ev' in kwargs) and ('Ec' in kwargs):
            self.Ev = kwargs['Ev']
            self.Ec = kwargs['Ec']
        
        else:
            self.Ev = [0 for i in range(nb_Ag)] # variable region binding energy for each antigen
            self.Ec = 0                         # constant region binding energy

            selected_Ag = np.random.randint(nb_Ag)
            for i in range(nb_Ag):
                if selected_Ag!=i: self.Ev[i] = o - np.exp(np.random.normal(mu, sigma))

        if 'nb_mut' in kwargs: self.nb_mut = kwargs['nb_mut']
        else:                  self.nb_mut = 0
        
        if 'last_bound' in kwargs: self.last_bound = kwargs['last_bound']
        else:                      self.last_bound = np.random.multinomial(self.nb, pvals = [1/float(nb_Ag)] * nb_Ag)

    @classmethod
    def clone(cls, b):
        return cls(1, Ev = [k for k in b.Ev], Ec = b.Ec, nb_mut = b.nb_mut, last_bound = b.last_bound) # Return a new copy of the input BCell
    
    def bind_to(self, Ag):
        """ Return binding energy with input antigen. """
        return self.Ec

    def divide(self):
        """ Run one round of division. """
        self.nb *= 2
    
    def mutate_CDR(self):
        """ Change in energy due to affinity-affecting CDR mutation. """
        #selected_Ag = np.random.randint(nb_Ag)
        #for i in range(nb_Ag):
        #    if selected_Ag!=i: self.Ev[i]  = np.random.rand()
        #    #else:              self.Ev[i] += o - np.exp(np.random.normal(mu, sigma))
        #    else:              self.Ev[i] += gevrnd.rvs()
        #self.Ec += o - np.exp(np.random.normal(mu, sigma))
        self.Ec     += gevrnd.rvs()
        self.nb_mut += 1

    def mutate_FR(self):
        """ Change in energy due to affinity-affecting framework (FR) mutation. """
        pass

    def shm(self):
        """
        Run somatic hypermutation and return self + new B cell clones.
        """
        
        # get number of cells that mutate
        new_clones = []
        n_mut      = np.random.binomial(self.nb, p_mut)
        self.nb   -= n_mut
            
        # get number of CDR vs framework (FR) mutations
        n_CDR = np.random.binomial(n_mut, p_CDR)
        n_FR  = n_mut - n_CDR
            
        # process CDR mutations
        n_die, n_silent, n_affect  = np.random.multinomial(n_CDR, pvals = [p_CDR_lethal, p_CDR_silent, p_CDR_affect])
        self.nb                   += n_silent
        for i in range(n_affect):
            b = BCell.clone(self)
            b.mutate_CDR()
            new_clones.append(b)
        
        # process FR mutations
        n_die, n_silent, n_affect  = np.random.multinomial(n_FR, pvals = [p_FR_lethal, p_FR_silent, p_FR_affect])
        self.nb                   += n_silent
        for i in range(n_affect):
            b = BCell.clone(self)
            b.mutate_FR()
            new_clones.append(b)

        # return the result
        if (self.nb>0): new_clones.append(self)
        return new_clones


###### Main functions ######


def usage():
    print("")


def main(verbose=False):
    """
    Simulate the affinity maturation process in a single germinal center (GC)
    and save the results to a CSV file.
    """
    
    # Run multiple trials and save all data to file
    
    nb_trial = 100
    start    = timer()
    
    fmem = open('output-memory.csv', 'w')
    ftot = open('output-total.csv',  'w')
    
    fmem.write('trial,exit cycle,number,mutations,Ec,'+(','.join(['Ev'+str(i) for i in range(nb_Ag)]))+'\n')
    ftot.write('trial,cycle,number recycled,number exit\n')
    
    for t in range(nb_trial):
    
        print_update(t, nb_trial)   # status check

        # INITIALIZATION - DEFINE DATA STRUCTURES

        memory_cells = []
        nb_recycled  = []
        nb_exit      = []

        
        # CYCLES 1 + 2 - CREATE FOUNDERS AND REPLICATE WITHOUT MUTATION
        
        nb_founders = 3                                     # number of founder B cells for a GC
        B_cells     = [BCell() for i in range(nb_founders)]
        
        # Update data
        nb_recycled.append(nb_founders)                     # all founders are recycled
        nb_exit.append(0)                                   # no founders exit the GC
        nb_recycled.append(np.sum([b.nb for b in B_cells])) # all founders replicate and are recycled
        nb_exit.append(0)                                   # no founders exit
        

        # AFFINITY MATURATION
        
        GC_size_max  = np.sum([b.nb for b in B_cells])  # maximum number of cells in the GC (= initial population size)
        nb_cycle_max = 250                              # maximum number of GC cycles
        cycle_number = 2
        
        for cycle_number in range(2, nb_cycle_max):
        
            B_cells, exit_cells = run_GC_cycle(B_cells)
            GC_size             = np.sum([b.nb for b in B_cells])       # total number of cells in the GC
            
            if (cycle_number==nb_cycle_max-1) or (GC_size>GC_size_max): # at the end, all B cells exit the GC
                exit_cells += B_cells
            else: exit_cells = []
            
            memory_cells.append(exit_cells)
            nb_recycled.append(np.sum([b.nb for b in B_cells]   ))
            nb_exit.append(    np.sum([b.nb for b in exit_cells]))

            if (nb_recycled[-1]==0) or (GC_size>GC_size_max): break
        

        # SAVE OUTPUT

        for i in range(len(memory_cells)):
            for b in memory_cells[i]:
                fmem.write('%d,%d,%d,%d,%lf' % (t, i+2, b.nb, b.nb_mut, b.Ec))
                for j in range(nb_Ag): fmem.write(',%lf' % b.Ev[j])
                fmem.write('\n')
        fmem.flush()

        for i in range(len(nb_recycled)): ftot.write('%d,%d,%d,%d\n' % (t, i+1, nb_recycled[i],nb_exit[i]))
        ftot.flush()

    # End and output total time
    
    fmem.close()
    ftot.close()
    
    end = timer()
    print('\nTotal time: %lfs, average per cycle %lfs' % ((end - start),(end - start)/float(nb_trial)))


def print_update(current, end, bar_length=20):
    """ Print an update of the simulation status. h/t Aravind Voggu on StackOverflow. """
    
    percent = float(current) / end
    dash    = ''.join(['-' for k in range(int(round(percent * bar_length)-1))]) + '>'
    space   = ''.join([' ' for k in range(bar_length - len(dash))])

    sys.stdout.write("\rSimulating: [{0}] {1}%".format(dash + space, int(round(percent * 100))))
    sys.stdout.flush()


def run_dark_zone(B_cells):
    """ B cells proliferate and undergo SHM in the dark zone. """
    
    n_rounds = 2
    for i in range(n_rounds):
        new_cells = []
        for b in B_cells:
            b.divide()
            new_cells += b.shm()
        B_cells = new_cells


def run_binding_selection(B_cells):
    """ Select B cells for binding to antigen. """
    
    for b in B_cells:
        
        b.last_bound = np.random.multinomial(b.nb, pvals = [1./float(nb_Ag)] * nb_Ag)
        
        for i in range(nb_Ag):
            
            # compute binding energy and chance of death ( = 1 - chance of survival )
            Ag_bound      = np.exp(energy_scale * b.bind_to(i))
            factor        = conc * Ag_bound
            langmuir_conj = 1. / (1. + factor)
            
            # remove dead cells and update binding details
            n_die            = np.random.binomial(b.last_bound[i], langmuir_conj)
            b.nb            -= n_die
            b.last_bound[i] -= n_die


def run_help_selection(B_cells):
    """ Select B cells to receive T cell help. """
    
    # get binding energies
    binding_energy     = [[b.bind_to(i) for i in range(nb_Ag)] for b in B_cells]
    binding_energy_tot = []
    for i in range(len(B_cells)):
        for j in range(nb_Ag): binding_energy_tot += [binding_energy[i][j]] * B_cells[i].last_bound[j]
    
    # cells in the top (help_cutoff) fraction of binders survive
    if len(binding_energy_tot)>0:
        cut_idx       = np.max([0, int(np.floor(help_cutoff * len(binding_energy_tot)))-1])
        energy_cutoff = np.array(binding_energy_tot)[np.argsort(binding_energy_tot)][::-1][cut_idx]
        n_die_tie     = len(binding_energy_tot) - cut_idx - np.sum(binding_energy_tot < energy_cutoff)

        # kill all B cells below threshold
        for i in np.random.permutation(len(B_cells)):
            for j in np.random.permutation(nb_Ag):
                energy = binding_energy[i][j]
                if energy < energy_cutoff:
                    B_cells[i].nb            -= B_cells[i].last_bound[j]
                    B_cells[i].last_bound[j]  = 0
                elif (energy == energy_cutoff) and (n_die_tie > 0):
                    if B_cells[i].last_bound[j] < n_die_tie:
                        B_cells[i].nb            -= B_cells[i].last_bound[j]
                        n_die_tie                -= B_cells[i].last_bound[j]
                        B_cells[i].last_bound[j]  = 0
                    else:
                        B_cells[i].nb            -= n_die_tie
                        B_cells[i].last_bound[j] -= n_die_tie
                        n_die_tie                 = 0


def run_recycle(B_cells):
    """ Randomly select B cells to be recycled back into the GC or to exit. """

    new_cells  = []                                 # cells that will remain in the GC
    exit_cells = []                                 # cells that will exit the GC
    n_tot      = np.sum([b.nb for b in B_cells])    # total number of cells currently in the GC
    n_exit     = int(np.floor(p_exit * n_tot))      # number of cells that will exit the GC
    b_exit     = []                                 # index of cells that will exit the GC

    if (n_tot > 0) and (n_exit > 0):
        b_exit = np.random.choice(n_tot, n_exit)

    idx = 0
    for b in B_cells:
    
        # find which cells exit the GC
        n_exit  = np.sum((idx <= b_exit) * (b_exit < idx + b.nb))
        idx    += b.nb
        b.nb   -= n_exit
        
        # add remainder to recycled cells
        if (b.nb>0):
            new_cells.append(b)
    
        # record exit cells
        if (n_exit>0):
            exit_cells.append(deepcopy(b))
            exit_cells[-1].nb = n_exit

    return new_cells, exit_cells

    # STOCHASTIC RECYCLING
#    for b in B_cells:
#        
#        n_exit  = np.random.binomial(b.nb, p_exit)
#        b.nb   -= n_exit
#        
#        if (b.nb>0):
#            new_cells.append(b)
#        
#        if (n_exit>0):
#            exit_cells.append(deepcopy(b))
#            exit_cells[-1].nb = n_exit
#
#    return new_cells, exit_cells


def run_GC_cycle(B_cells):
    """ Run one cycle of the GC reaction. """

    run_dark_zone(B_cells)          # DARK  ZONE - two rounds of division + SHM
    run_binding_selection(B_cells)  # LIGHT ZONE - selection for binding to Ag
    run_help_selection(B_cells)     # LIGHT ZONE - selection to receive T cell help
    return run_recycle(B_cells)     # RECYCLE    - randomly pick exiting cells from the surviving B cells


################################################

def test_dark_zone():
    """ Test the run_dark_zone function. """
    
    n_clones    = []
    n_cells_max = []
    n_tests     = 10000
    
    for i in range(n_tests):
        test_cells = [BCell(nb = 1000)]
        run_dark_zone(test_cells)

        temp_clones = 0
        temp_max    = 0
        for b in test_cells:
            if b.nb>0:        temp_clones += 1
            if b.nb>temp_max: temp_max=b.nb
        n_clones.append(temp_clones)
        n_cells_max.append(temp_max)

    



def run_tests():
    """ Run diagnostic tests to make sure that the code is functioning as expected. """
    pass



if __name__ == '__main__': main()

