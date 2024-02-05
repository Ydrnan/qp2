#!/bin/python3
from classes import *

def get_lines_match(f_path, match):
    acc = []
    with open(f_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if match in line:
                acc.append(line.strip())

    return acc

def go_to_match(buff,i,match):
    while i < len(buff):
        if not match in buff[i]:
            i += 1
            continue
        else:
            break

    if i >= len(buff):
        print("Reaches end of file.")
        sys.exit()

    return i

def try_go_to_match(buff,i,match):
    while i < len(buff):
        if not match in buff[i]:
            i += 1
            continue
        else:
            break

    return i

if __name__ == "__main__":
    import sys
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    files = sys.argv[1:]

    l_cipsi = []
    for fname in files:
        cipsi = Cipsi([])
        with open(fname,'r') as f:
            buff = f.readlines()
            i = 0

            n_states = int(get_lines_match(fname,"Number of states:")[0].split()[-1])
            n_eta = int(get_lines_match(fname,"Number of eta values:")[0].split()[-1])

            first = True
            while i < len(buff):

                i = go_to_match(buff,i,"Number of eta values:")

                # Jump to the values
                i += 8
                buff_eta = buff[i:i+n_eta]
                a = []
                for line in buff_eta:
                    for val in line.split():
                        a.append(val)
                a = np.asarray(a, dtype=float, order='F')
                a = a.reshape((n_eta,2*n_states+1))

                # N_det
                i = go_to_match(buff,i,"Summary at N_det =")    
                n_det = int(buff[i].split()[-1])

                # E                    
                i = go_to_match(buff,i,"# E ")
                l_str_E = buff[i].split()[2:n_states+2]
                l_E = [float(e) for e in l_str_E]

                # PT2
                i += 3
                l_str_pt2 = buff[i].split()[2:2*n_states+2]
                l_pt2 = [float(e) for e in l_str_pt2[::2]]

                l_states = []
                l_eta = a[:,0]
                for j,E,pt2 in zip(range(1, 2*n_states+1, 2),l_E,l_pt2):
                    l_energy = []
                    l_re = a[:,j]
                    l_im = a[:,j+1]
                    for eta, re, im in zip(l_eta,l_re,l_im):
                        l_energy.append(Energy(float(eta),float(re),float(im)))
 
                    energies = Energies(n_eta,l_energy)
                    state = (j+1)//2

                    l_states.append(State(n_det,E,pt2,energies,state))

                step = Step(n_states, l_states)

                cipsi.append(step)

                i = try_go_to_match(buff,i,"Number of eta values:")

        l_cipsi.append(cipsi)

    data = []
    for cispi in l_cipsi:
        for step in cipsi.steps:
            for state in step.states:
                print(state.n_det,state.E,state.pt2,state.state)
                for e in state.energies.energies:
                    print(e.eta, e.re, e.im)
                    
                    data.append([state.state, e.eta, e.re, e.im, state.n_det, state.E, state.pt2])

    df = pd.DataFrame(data, columns=['state','eta','Re(E_cap)','Im(E_cap)','N_det','E','PT2'])
    
    print(df)
    N_det = df['N_det'].unique().tolist()
    states = df['state'].unique().tolist()
    val_eta = df['eta'].unique().tolist()

    print(N_det, states, val_eta)

    last_pt_df = df[(df['N_det'] == 7473) & (df['state'] == 2)]    

    print(last_pt_df)

    last_pt_df.plot(kind = 'scatter', x = 'eta', y = 'Re(E_cap)')
    plt.show()

