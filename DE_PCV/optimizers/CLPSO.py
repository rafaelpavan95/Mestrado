
import numpy as np

class Comprehensive_Learning_PSO():
    
    '''
    Creates the Comprehensive Learning Particle Swarm Optimization Algorithm.
    '''
    ############################################
    
    def __init__(self, refresh, c1, c2, zeta, gamma, w_max, w_min, max_iter, n_particles, alfa, beta ,v_amp, v_clamp, initial=None, initi = False):
        
        '''
        Initializes the hyperparameters for the metaheuristic.
        '''
        self.c1=c1 # learning factor for cognitive component
        self.c2=c2
        self.refresh_gap = refresh
        self.max_iter = max_iter # algorithm maximum number of iterations
        self.n_particles = n_particles # algorithm number of particles
        self.alfa = alfa # penalization factor for costs
        self.beta = beta # penalization factor for demand mismatch
        self.v_amp = v_amp # amplitude atenuation for velocity
        self.v_clamp = v_clamp # clamping atenuation for velocity 
        self.initial = initial
        self.zeta = zeta
        self.gamma = gamma
        self.flag = initi
        self.wmax = w_max
        self.wmin = w_min

    ############################################
    def init_pop(self, system):

        '''
        Initializes the swarm (randomly).
        '''

        population = np.zeros((self.n_particles,len(system.pg_min)))
        
        for gen in range(len(system.pg_min)):
            
            population[:,gen] = np.random.uniform(system.pg_min[gen],system.pg_max[gen],size=([self.n_particles]))
    
        if self.flag == True:
            
            population[:,:]=self.initial
            
        return population
    
    ############################################
    def init_vel(self,system):

        '''
        Initializes the velocities (randomly).
        '''

        vel = np.zeros((self.n_particles,len(system.pg_min)))
        
        for gen in range(len(system.pg_min)):
            
            vel[:,gen] = np.random.uniform(-system.pg_max[gen],system.pg_max[gen],size=([self.n_particles]))
        
        return vel

    ############################################

    
    ############################################
    def objective_function(self, system, population):

        evaluated = np.zeros(self.n_particles)
        
        for particle in range(self.n_particles):
            
            evaluated[particle] = system.costs(population[particle,:])
        
        return evaluated
    
    ############################################
    def inequality_constraints(self, system, population):
        
        for gen in range(len(system.pg_min)):
            
            for part in range(self.n_particles):

                if population[part,gen] < system.pg_min[gen]:
                    population[part,gen] = system.pg_min[gen]
                if population[part,gen] > system.pg_max[gen]:
                    population[part,gen] = system.pg_max[gen]

        return population
           
    ############################################
    def equality_constraints(self, system, population):
        
        generation = np.sum(population,axis=1)
        
        balance = (generation-system.demanda)      
        
        # balance[balance>0] = 0
        
        balance = np.abs(balance)
        
        return balance
    
    ############################################
    def fitness_function(self, system, costs, balance):
            
        fitness = self.alfa*costs + self.beta*balance
        
        return fitness
    
    ############################################
    def find_gbest(self,fitness,population):
        
        position_gbest = np.argmin(fitness)
        
        gbest = np.copy(population[position_gbest,:])
        
        return gbest, np.min(fitness)
    
    ############################################
    def update_memory(self, system, fitness, population, memory):
        
         
        ev_mem = self.objective_function(system, memory)
        eq_mem = self.equality_constraints(system, memory)
        fit_mem = self.fitness_function(system, ev_mem, eq_mem)
            
        for i in range(len(fit_mem)):
                
            if fit_mem[i] > fitness[i]:
                    
                memory[i,:] = np.copy(population[i,:])
        
        return memory
    
    ############################################
    def update_inertia(self, system, iteration):
  
        w = self.wmax - (self.wmax-self.wmin)*iteration/self.max_iter
        
        return w
    
    ############################################
    def pbestfm(self, proba, memoria, system, flag, pbestfcopy):
        
        ev_mem = self.objective_function(system, memoria)
        eq_mem = self.equality_constraints(system, memoria)
        fit_mem = self.fitness_function(system, ev_mem, eq_mem)
 
    
        pbestf = np.copy(memoria)
        
        for i in range(self.n_particles):
            if flag[i] >=self.refresh_gap:
                
                flag[i]=0
                for d in range(len(system.pg_min)):

                    aleatorio = np.random.rand()


                    if aleatorio<proba[i]:

                        p1 = int(np.random.rand()*self.n_particles)
                        p2 = int(np.random.rand()*self.n_particles)

                        while((p2==i) or (p1==i) or (p1==p2)):

                            p1 = int(np.random.rand()*self.n_particles)
                            p2 = int(np.random.rand()*self.n_particles)

                        if fit_mem[p1]<fit_mem[p2]:

                            pbestf[i,d] = memoria[p1,d]
                            
                        else:

                            pbestf[i,d] = memoria[p2,d]
                          


                    else:

                        pbestf[i,d] = memoria[i,d]
                        
            else:
                
                pbestf[i,:] = pbestfcopy[i,:]


        for i in range(self.n_particles):

            summ = 0
            
            for d in range(len(system.pg_min)):
            
                if pbestf[i,d] == memoria[i,d]:
                    
                    summ = summ+1
                    
            if summ == len(system.pg_min):
                
                pt = int(np.random.rand()*self.n_particles)
                d = int(np.random.rand()*len(system.pg_min))
                
                pbestf[i,d] = memoria[pt,d]
                           

        return pbestf
    
    ############################################
    def update_velocity(self, previous_v, w, population, pbestf, gbest,system):

        r1 = np.random.rand(self.n_particles,len(system.pg_max))
        r2 = np.random.rand(self.n_particles,len(system.pg_max))
        
        new_v = previous_v*w + self.c1*r1*(pbestf-population) + self.c2*r2*(gbest-population)
        
        return new_v

    ############################################
    def velocity_clamping(self, system, v):
        
        delta = np.abs(system.pg_max-system.pg_min)
        
        for j in range(self.n_particles):
            
            for k in range(len(system.pg_min)):
                
                if v[j,k] < -delta[k]*self.v_clamp:
                    
                    v[j,k] = -delta[k]*self.v_clamp
                    
                
                if v[j,k] > delta[k]*self.v_clamp:
                    
                    v[j,k] = delta[k]*self.v_clamp
        
        return v

    ############################################
    def update_position(self,new_v,population,system):
        
           
        new_position = np.copy(population) + new_v
        
        
        return new_position
    
    ############################################

    def adjust_dem(self,a, system):

        DIF = system.demanda-np.sum(a)

        parcela = DIF/len(system.pg_min)

        a = a + parcela

        return a

    ############################################
    def adjust_viol(self,a, system):

        flag = np.ones(len(system.pg_min))

        for gen in range(len(system.pg_min)):

            if a[gen] < system.pg_min[gen]:
                a[gen] = system.pg_min[gen]

            elif a[gen] > system.pg_max[gen]:
                a[gen] = system.pg_max[gen]
            else:
                flag[gen]=0

        return a, flag

    ############################################
    def adjust(self,a,flag,system, fator):

        while (np.abs((np.sum(a)-system.demanda))>fator and np.sum(flag)!=0):

            a, flag = self.adjust_viol(a,system)

            a = self.adjust_dem(a,system)

            a, flag = self.adjust_viol(a,system)

        return a

    ############################################
 
    def optimize(self, system, verbose=True):
        
        pop = self.init_pop(system)

        previous_vel = self.init_vel(system)*self.v_amp
        
        self.list_gbest = []
        self.list_gbest_costs = []
        self.list_gbest_pen = []
        self.list_gbest_fitness = []
        self.list_costs_evaluation = []
        self.list_generation_pen = []
        self.list_fitness = []
        self.list_inertia = []
        self.global_bests = []
        self.global_bests_val = []

        swarm = []
        
        refresh_gap = np.ones(self.n_particles)*self.refresh_gap
        
        pbestfc = np.copy(pop)*0
        
        for iterations in range(self.max_iter):

            swarm.append(pop)
            
            costs = self.objective_function(system, np.copy(pop))

            self.list_costs_evaluation.append(costs)
            
            generation_balance = self.equality_constraints(system, np.copy(pop))

            self.list_generation_pen.append(generation_balance)
            
            fit = self.fitness_function(system, costs, generation_balance)
            
            self.list_fitness.append(fit)
            
            gbest, gbest_val = self.find_gbest(fit,np.copy(pop))

            flag = np.ones(len(system.pg_min))
            
            if iterations>self.max_iter*0.8:
                
                gbest = self.adjust(gbest,flag,system,  1e-3)

  
            if iterations == 0:
                
                self.global_bests.append(gbest)

                self.global_bests_val.append(gbest_val)
                

            else:
                
                if self.global_bests_val[-1] > gbest_val:
                    # from optimizers.SQP import Sequential_Quadratic_Programming
                    # sqp = Sequential_Quadratic_Programming()

                    # r = sqp.optimize(system, gbest.copy(), verbose=False)

                    # if r['gbest_cost'] < gbest_val:

                    #     gbest_val=r['gbest_cost']
                    #     gbest = r['gbest']

                    self.global_bests_val.append(gbest_val)
                    self.global_bests.append(gbest)
          
            
            self.list_gbest.append(self.global_bests[-1])

            self.list_gbest_pen.append(self.equality_constraints(system, self.global_bests[-1].reshape(1,len(system.pg_min))))
            
            self.list_gbest_costs.append(system.costs(self.global_bests[-1]))
            
            self.list_gbest_fitness.append(self.fitness_function(system, self.list_gbest_costs[-1], self.list_gbest_pen[-1])[0])
            
            if iterations == 0:
            
                pbest = np.copy(pop)
            
            else:

                pbest = self.update_memory(system, fit, np.copy(pop), np.copy(pbest))
                
            
            w = self.update_inertia(system, iterations)
            
            self.list_inertia.append(w)
            
            p = np.arange(start=0,stop=self.n_particles+1,step=1)
            proba = self.zeta + self.gamma*((np.exp((10*(p-1))/(30-1))-1)/(np.exp(10)-1))
            proba = proba[1:]
            pbestf = self.pbestfm(proba, np.copy(pbest), system, refresh_gap,pbestfc)
            pbestfc = np.copy(pbestf)
  
            if iterations>1:
                
                for i in range(self.n_particles):
                    
                    if self.list_fitness[-1][i]>=self.list_fitness[-2][i]:
                        refresh_gap[i] = refresh_gap[i]+1
                

            new_vel = self.update_velocity(previous_vel, w, np.copy(pop), pbestf, gbest, system)
            
            new_vel = self.velocity_clamping(system, new_vel)
            
            previous_vel = np.copy(new_vel)
            
            new_position = self.update_position(new_vel,np.copy(pop),system)
            
            new_position = self.inequality_constraints(system, np.copy(new_position))
            
            pop = np.copy(new_position)
            
            if verbose==True:
            
                print(f'Iteration: {iterations}')
                print(f'Inertia Weight: {w}')
                print(f'Gbest Costs [$]: {self.list_gbest_costs[-1]}')
                print(f'Gbest Penalization: {self.list_gbest_pen[-1][0]}')
                print(f'Gbest Fitness: {self.list_gbest_fitness[-1]}')
                print(f'Gbest Active Power Generation [MW]: {np.sum(self.list_gbest[-1])}')
                print('_ _ _ _ _ _ _ _ _ _ _ \n')
        
        return {'gbest_list': self.list_gbest,
                'gbest_costs_list': self.list_gbest_costs,
                'gbest_pen_list': self.list_gbest_pen,
                'gbest_list_fitness': self.list_gbest_fitness}




