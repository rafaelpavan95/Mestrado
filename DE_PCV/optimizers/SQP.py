import numpy as np
from gekko import GEKKO

class Sequential_Quadratic_Programming():
    
    def obj(self,pg,a_k,b_k,c_k,e_k,f_k,pg_min,pg_max,alfa,model):
    
        soma = 0
        for i in range(len(pg_min)):

            val = e_k[i]*model.sin(f_k[i]*(pg_min[i]-pg[i]))

            soma = pg[i]*pg[i]*a_k[i] + pg[i]*b_k[i] + c_k[i] + alfa[i] + soma

        return soma
    
        
    def constraint(self,pg):

        soma = 0

        for i in range(len(pg)):

            soma = pg[i] + soma

        return soma
        
    
    def optimize(self, system, initial, verbose=True):
        
        gbest = initial
        model = GEKKO(remote=False)

        pg=[]
        pg_min = []
        pg_max = []
        a_k = []
        b_k = []
        c_k = []
        e_k = []
        f_k = []
        alfa = []
        demand = system.demanda
        
        for i in range(len(system.pg_min)):
            
                            
            pg.append(model.Var(gbest[i],system.pg_min[i],system.pg_max[i]))

            pg_min.append(system.pg_min[i])

            pg_max.append(system.pg_max[i])

            a_k.append(system.a_k[i])

            b_k.append(system.b_k[i])

            c_k.append(system.c_k[i])

            e_k.append(system.e_k[i])

            f_k.append(system.f_k[i])
            
            alfa.append(model.Var(np.abs(e_k[i]*np.sin(f_k[i]*(system.pg_min[i]-gbest[i]))),0,1e15))

        model.Equation(self.constraint(pg)==demand)
        
        ind=0
        
        for val in alfa:
            
            model.Equation(val>=e_k[ind]*model.sin(f_k[ind]*(system.pg_min[ind]-pg[ind])))

            model.Equation(val >= -e_k[ind]*model.sin(f_k[ind]*(system.pg_min[ind]-pg[ind])))

            ind=ind+1
        
        model.options.CSV_WRITE=2

        model.Obj(self.obj(pg,a_k,b_k,c_k,e_k,f_k,pg_min,pg_max,alfa,model))
        
        model.options.SOLVER=1

        model.solver_options = ['constraint_convergence_tolerance 1.0e-3']

        model.solve(disp=verbose)
        
        
        
        return {'gbest': np.asarray(pg).ravel(), 'gbest_cost':system.costs(np.asarray(pg).ravel())}
        
        
