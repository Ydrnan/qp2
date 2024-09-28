#/bin/python3
import sys

class Energy:
    def __init__(self,eta,E_re,E_im,C_re,C_im):
        if type(eta) != float or type(E_re) != float or type(E_im) != float \
        or type(C_re) != float or type(C_im) != float:
            print("Invalid type when creating Energy class object.")
            sys.exit()
            
        self.eta = eta
        self.re = E_re
        self.im = E_im
        self.c_re = C_re
        self.c_im = C_im

class Energies:
    def __init__(self,n_eta, l_energy):
        if type(l_energy) != list:
            print("Invalid type when creating Energies class object.")
            sys.exit()

        for e in l_energy:
            if type(e) != Energy:
                print("Invalid type, Energies must take a list of Energy objects.")
                sys.exit()

        if len(l_energy) != n_eta:
            print("Wrong number of eta values")
            sys.exit()

        self.n_eta = n_eta
        self.energies = l_energy

class State:
    def __init__(self,n_det,E,pt2,energies,state):
        if type(n_det) != int:
            print("Invalid type for n_det when creating State class object.")
            sys.exit()
        if type(E) != float:
            print("Invalid type for E when creating State class object.")
            sys.exit()
        if type(pt2) != float:
            print("Invalid type for pt2 when creating State class object.")
            sys.exit()
        if type(state) != int:
            print("Invalid type for state when creating State class object.")
            sys.exit()

        self.n_det = n_det
        self.E = E
        self.pt2 = pt2
        self.energies = energies
        self.state = state

class Step:
    def __init__(self, n_states, states):
        if type(states) != list:
            print("Invalid type when creating Step class object.")
            sys.exit()

        if len(states) != n_states:
            print("Wrong number of states")
            sys.exit() 

        self.n_states = n_states
        self.states = states

class Cipsi: 
    def __init__(self,steps):
        if type(steps) != list:
            print("Invalid type, Cipsi class objects are creating with a list of Step.")
            sys.exit()

        self.steps = steps                

    def append(self,step):
        if type(step) != Step:
            print("Invalid type, append methods in Cipsi must take an Step class object.")
            sys.exit()

        self.steps.append(step)


