import numpy as np

class Particle_Swarm_Optimization():

    def __init__(self, c1_, c2_, w_max, w_min, inertia_method, max_iter_, n_particles_, alfa_, beta_ ,v_amp_, v_clamp_, initial_, initi = True, constriction_factor = False, inertia_weight=True):
        
        self.c1=c1_ # learning factor forg cognitive component
        self.c2=c2_ # learning factor for social component
        self.max_iter = max_iter_ # algorithm maximum number of iterations
        self.n_particles = n_particles_ # algorithm number of particles
        self.alfa = alfa_ # penalization factor for costs
        self.beta = beta_ # penalization factor for demand mismatch
        self.v_amp = v_amp_ # amplitude atenuation for velocity
        self.v_clamp = v_clamp_ # clamping atenuation for velocity 
        self.initial = initial_
        
        self.flag = initi
        
        if constriction_factor == True:
            
            pass
        
            if self.c1+self.c2 <= 4:
                print('Error: c1+c2 must be >= 4')
                
            
        else:
            
            self.wmax = w_max
            self.wmin = w_min
            self.inertia = str(inertia_method)

    def init_pop(self, system):
        
        population = np.zeros((self.n_particles,len(system.pg_min)))
        
        for gen in range(len(system.pg_min)):
            
            population[:,gen] = np.random.uniform(system.pg_min[gen],system.pg_max[gen],size=([self.n_particles]))
    
        if self.flag == True:
            population[0,:]=self.initial
            
        return population

    def init_vel(self,system):
        
        vel = np.zeros((self.n_particles,len(system.pg_min)))
        
        for gen in range(len(system.pg_min)):
            
            vel[:,gen] = np.random.uniform(-system.pg_max[gen],system.pg_max[gen],size=([self.n_particles]))
        
        return vel

    def objective_function(self, system, population):
        
        evaluated = np.zeros(self.n_particles)
        
        for particle in range(self.n_particles):
            
            evaluated[particle] = system.costs(population[particle,:])
        
        return evaluated

    def inequality_constraints(self, system, population):
        
        for gen in range(len(system.pg_min)):
            
            for part in range(self.n_particles):

                if population[part,gen] < system.pg_min[gen]:
                    population[part,gen] = system.pg_min[gen]
                if population[part,gen] > system.pg_max[gen]:
                    population[part,gen] = system.pg_max[gen]

        return population

    def equality_constraints(self, system, population):
        
        generation = np.sum(population,axis=1)
        
        balance = (generation-system.demanda)      
        
        balance = np.abs(balance)
        
        return balance

    def fitness_function(self, system, costs, balance):
            
        fitness = self.alfa*costs + self.beta*balance
            
        return fitness

    def find_gbest(self,fitness,population):
        
        position_gbest = np.argmin(fitness)
        
        gbest = np.copy(population[position_gbest,:])
        
        return gbest, np.min(fitness)

    def update_memory(self, system, fitness, population, memory):
        
         
        ev_mem = self.objective_function(system, memory)
        eq_mem = self.equality_constraints(system, memory)
        fit_mem = self.fitness_function(system, ev_mem, eq_mem)
            
        for i in range(len(fit_mem)):
                
            if fit_mem[i] > fitness[i]:
                    
                memory[i,:] = np.copy(population[i,:])
        
        return memory

    def update_inertia(self, system, iteration):
        
        if self.inertia == 'linear':
            
            w = self.wmax - (self.wmax-self.wmin)*iteration/self.max_iter
            
        if self.inertia == '':
        
            w = 1
        
        return w

    def update_velocity(self, previous_v, w, population, memory, gbest, system):
        
        r1 = np.random.rand(self.n_particles,len(system.pg_max))
        r2 = np.random.rand(self.n_particles,len(system.pg_max))
        
        
        new_v = previous_v*w + self.c2*r2*(gbest-population) + self.c1*r1*(memory-population)
        
        return new_v

    def velocity_clamping(self, system, v):
        
        delta = np.abs(system.pg_max-system.pg_min)
        
        for j in range(self.n_particles):
            
            for k in range(len(system.pg_min)):
                
                if v[j,k] < -delta[k]*self.v_clamp:
                    
                    v[j,k] = -delta[k]*self.v_clamp
                    
                
                if v[j,k] > delta[k]*self.v_clamp:
                    
                    v[j,k] = delta[k]*self.v_clamp
        
        return v

    def update_position(self,new_v,population, const_factor):
        
        new_position = np.copy(population) + new_v*1
        
        return new_position
    
    def adjust_dem(self,a, system):

        DIF = system.demanda-np.sum(a)

        parcela = DIF/len(system.pg_min)

        a = a + parcela

        return a


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


    def adjust(self,a,flag,system, fator):

        while (np.abs((np.sum(a)-system.demanda))>fator and np.sum(flag)!=0):

            a, flag = self.adjust_viol(a,system)

            a = self.adjust_dem(a,system)

            a, flag = self.adjust_viol(a,system)

        return a


    def optimize(self, system, verbose=True):
        
        pop = self.init_pop(system)
        
        previous_vel = self.init_vel(system)*self.v_amp
        
        self.list_gbest = []
        
        self.list_gbest_costs = []
        self.list_gbest_pen = []
        self.list_gbest_fitness = []
        
        self.list_pop = []
        self.list_costs_evaluation = []
        self.list_generation_pen = []
        self.list_fitness = []
        self.list_inertia = []
        swarm = []
        global_bests = []
        global_bests_val = []
        
        for iterations in range(self.max_iter):

            swarm.append(pop)
            
            self.list_pop.append(np.copy(pop))
            
            costs = self.objective_function(system, np.copy(pop))
            
            self.list_costs_evaluation.append(costs)
            
            generation_balance = self.equality_constraints(system, np.copy(pop))
            
            self.list_generation_pen.append(generation_balance)
            
            fit = self.fitness_function(system, costs, generation_balance)
            
            self.list_fitness.append(fit)
            
            gbest, gbest_val = self.find_gbest(fit,np.copy(pop))
            
            flag = np.ones(len(system.pg_min))
            
            if iterations>self.max_iter*0.8:
                gbest = self.adjust(gbest,flag,system, 1e-3)
      
            if iterations == 0:
                
                global_bests.append(gbest)
                global_bests_val.append(gbest_val)
                
            else:
                
                if global_bests_val[-1] > gbest_val:
                    
                    global_bests_val.append(gbest_val)
                    global_bests.append(gbest)
          
            self.list_gbest.append(global_bests[-1])
            
            self.list_gbest_pen.append(self.equality_constraints(system, global_bests[-1].reshape(1,len(system.pg_min))))
            
            self.list_gbest_costs.append(system.costs(global_bests[-1]))
            
            self.list_gbest_fitness.append(self.fitness_function(system, self.list_gbest_costs[-1],self.list_gbest_pen[-1]))
            
            if iterations == 0:
            
                pbest = np.copy(pop)
            
            else:

                pbest = self.update_memory(system, fit, np.copy(pop), np.copy(pbest))
                
            
            w = self.update_inertia(system, iterations)
            
            self.list_inertia.append(w)
            
            new_vel = self.update_velocity(previous_vel, w, np.copy(pop), np.copy(pbest), global_bests[-1], system)
            
            new_vel = self.velocity_clamping(system, new_vel)
            
            previous_vel = np.copy(new_vel)
            
            new_position = self.update_position(new_vel,np.copy(pop), 1)
            
            new_position = self.inequality_constraints(system, np.copy(new_position))
            
            pop = np.copy(new_position)
            
            if verbose==True:
            
                print(f'Iteration: {iterations}')
                print(f'Inertia Weight: {w}')
                print(f'Gbest Costs [$]: {self.list_gbest_costs[-1]}')
                print(f'Gbest Penalization: {self.list_gbest_pen[-1][0]}')
                print(f'Gbest Fitness: {self.list_gbest_fitness[-1][0]}')
                print(f'Gbest Active Power Generation [MW]: {np.sum(self.list_gbest[-1])}')
                print(f'Penalização: {self.beta}')
                print('_ _ _ _ _ _ _ _ _ _ _ \n')
        
        return self.list_gbest[-1]


        
