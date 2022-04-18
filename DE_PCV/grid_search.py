import numpy as np
import pandas as pd
from Systems import Economic_Dispatch
from Hybrid_CLPSO_SQP import Comprehensive_Learning_PSO
import time

def exaustive_search(n_gen, iterations):

    system_tests = {'refresh_gap': [0,1,2,3], 'c1': [1,2,3], 'zeta': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]}

    table = pd.DataFrame({'refresh_gap':[], 'c':[], 'probability':[], 'mean':[], 'min':[], 'std':[]}) 

    refresh_ = []

    c_ = []

    zeta_ = []

    mean_ = []

    min_ = []

    std_ = []

    prob_ = []

    values = []

    system = Economic_Dispatch(n_gen)

    for r_gap in system_tests['refresh_gap']:

        print('refresh gap: ', r_gap)

        for c in system_tests['c1']:
            
            print('learning_rate: ', c)

            for z in system_tests['zeta']:

                print('probability: ', z)    

                for i in range(iterations):

                    CLPSO = Comprehensive_Learning_PSO(refresh=r_gap, c1=c, c2=0, zeta=z, gamma=0, w_max=0.9, w_min=0.4, max_iter=1000, n_particles=100, alfa=1, beta=1e3 ,v_amp=0.1, v_clamp=0.1, initial=None, initi = False)
                    data = CLPSO.optimize(system, verbose=False)
                    values.append(data['gbest_list_fitness'][-1])
                    print(data['gbest_list_fitness'][-1])
                    print(i)
                    
                refresh_.append(r_gap)
                c_.append(c)
                prob_.append(z)
                mean_.append(np.mean(values))
                min_.append(np.min(values))                
                std_.append(np.std(values))



    table['refresh_gap'] = refresh_
    table['c'] = c_
    table['mean'] = mean_
    table['min'] = min_
    table['std'] = std_
    table['probability'] = prob_

    table.to_csv('tests_13_de_pcv.csv')

    return table



def automated_tests_clpso(system, d, max_ex, verb=False):

    '''
    input:
    system -> system to performe optimization process
    d -> dictionary with best parameters selected using gridsearch
    max_ex -> number of executions

    output:
    solution for each execution
    convergence curves for each execution
    '''

    clpso = Comprehensive_Learning_PSO(refresh=d['refresh'], c1=d['c1'], c2=0, zeta=d['probability'], gamma=0, w_max=0.9, w_min=0.4, max_iter=d['max_iterations'], n_particles = d['n_particles'], alfa=1, beta=50, v_amp=d['amp'], v_clamp=d['clamp'], initial=None, initi=False)

    solutions_clpso = {}

    convergence_curves = {}

    time_ = []

    for i in range(max_ex):
        start = time.time()
        results = clpso.optimize(system, verbose=verb)
        solutions_clpso[str(i)] = results['gbest_list'][-1]
        convergence_curves[str(i)] = results['gbest_list_fitness']
        end = time.time()
        time_.append(end-start)

    
    print('Average Time (s): ', np.mean(time_))
    print('Standard Deviation Time (s): ', np.std(time_))
    return pd.DataFrame(solutions_clpso), pd.DataFrame(convergence_curves)

