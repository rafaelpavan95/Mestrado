from resources.msc_rafael_pavan import *

def Particle_Swarm_Optimization(sep, zeta, psi, sigma, omega, neta,max_iter, n_particles,c1,c2,v_amp,valor_inicial,step, wmax,wmin,relatorio=True,inicial=True):
        
    enxame_fit = cria_enxame_fpo(sep,n_particles)

    if len(sep.bus) == 14:
        
        n_vgen = 4+1
        n_tap = 3
        n_gens = 4
        n_bshunt = 1
    
    if len(sep.bus) == 30:
        
        n_vgen = 5+1
        n_tap = 4
        n_bshunt = 2
        n_gens = 5
        
    
    if len(sep.bus) == 118:
        
        n_vgen = 53+1
        n_tap = 9
        n_bshunt = 14
        n_gens=53
        
        
    if len(sep.bus) == 300:
        
        n_vgen = 68+1
        n_tap = 62
        n_bshunt = 29
        n_gens=68
        
    w_max=wmax
    
    w_min=wmin
    
    
    j = []
    
    
    tempo = []
        
    obj = []
    
    pen_v = []
    
    pen_gq = []
    
    pen_tap = []
    
    pen_bsh = []
    
    pen_slack = []

    
    v_lim_superior = np.repeat(sep.bus['max_vm_pu'][0], n_vgen)
    
    v_lim_inferior = np.repeat(sep.bus['min_vm_pu'][0], n_vgen)
    
    tap_pos, tap_neutral, tap_step_percent,valores_taps = coleta_dados_trafo(sep,relatorio=False)
    
    tap_max = np.repeat(valores_taps[-1], len(tap_pos))
    
    tap_min = np.repeat(valores_taps[0], len(tap_pos))
    
    max_pot = sep.gen['max_p_mw'].values/100
    
    min_pot = sep.gen['min_p_mw'].values/100
    
    bsh,b=coleta_dados_bshunt(sep)

    bsh_max=[]
    
    bsh_min=[]
       

    for bs in bsh:
        bsh_max.append([np.max(bs)])
        bsh_min.append([np.min(bs)])


    maximo = np.expand_dims(np.concatenate((v_lim_superior, tap_max, bsh_max, max_pot), axis = None), 0)
    minimo = np.expand_dims(np.concatenate((v_lim_inferior, tap_min, bsh_min, min_pot), axis = None), 0)
     
    
    lim_sup = np.tile(maximo, (n_particles,1))
    lim_inf = np.tile(minimo, (n_particles,1))
    
    v_anterior = v_amp*cria_enxame_fpo(sep,n_particles)
    
    delta= np.abs(lim_sup-lim_inf)[0,:]
    

    if inicial == True:
        
        enxame_fit[i,:]=valor_inicial
            
    for i in range(0,max_iter):
        
            
        start = time.time()
        
        r1 = np.random.random_sample(size = (n_particles,enxame_fit.shape[1]))
        
        r2 = np.random.random_sample(size = (n_particles,enxame_fit.shape[1]))

        enxame_fit_d = np.copy(enxame_fit)
    
        for linha in range(n_particles):
          
            enxame_fit_d[linha][n_vgen:n_vgen+n_tap] = discreto_tap(enxame_fit[linha].copy(),n_tap,n_vgen,n_bshunt,sep)
            enxame_fit_d[linha][n_vgen+n_tap:n_vgen+n_tap+n_bshunt] = discreto_bshunt(enxame_fit[linha].copy(),n_tap,n_vgen,n_bshunt,sep)
            enxame_fit_d[linha][n_vgen+n_tap+n_bshunt:n_vgen+n_tap+n_bshunt+n_gens] = poz(n_vgen, n_tap, n_bshunt, n_gens, enxame_fit[linha].copy(), sep)[n_vgen+n_tap+n_bshunt:n_vgen+n_tap+n_bshunt+n_gens] 
            
        enxame_fit[:,-7:] = (fluxo_de_pot_fpo(enxame_fit_d,sep))[:,-7:]
        enxame_fit[:,-7:] = (fitness_fpo(enxame_fit,zeta,psi,sigma,omega,neta))[:,-7:]
        
        
        if i==0:
            
            best_particles = enxame_fit.copy()

            global_best = best_particles[np.argsort(best_particles[:, -1])][0,:].copy()
            
            global_matriz = np.tile(global_best, (n_particles,1))
        

        for t in range(0,n_particles):
                
            if (enxame_fit[t,-1] < best_particles[t,-1]):
        
                best_particles[t,:] = enxame_fit[t,:].copy()
            

        global_best = best_particles[np.argsort(best_particles[:, -1])][0,:].copy()

        global_matriz = np.tile(global_best, (n_particles,1))   
            
        enxame_fit_anterior = enxame_fit.copy()
            
        w_novo = w_max-(w_max-w_min)*(i+1)/max_iter 
                            
        v_novo = np.multiply(w_novo,v_anterior.copy()) + c1*np.multiply(r1,(best_particles.copy()-enxame_fit.copy())) + c2*np.multiply(r2,(global_matriz.copy()-enxame_fit.copy()))
               
        for l in range(v_novo.shape[0]):
        
            for c in range(len(delta)):
                
                if v_novo[l,c] > delta[c]*step:
                    
                    v_novo[l,c] = delta[c].copy()*step
                    
                    
                if v_novo[l,c] < -delta[c]*step:
                    
                    v_novo[l,c] = -delta[c].copy()*step
                    
        
                            
        enxame_fit_novo = enxame_fit_anterior  + v_novo        
    
        v_anterior = v_novo.copy()           
        
        enxame_estat = enxame_fit_novo[:,-7:]

        enxame_fit = np.concatenate(( np.clip(enxame_fit_novo[:,0:-7], a_min = lim_inf, a_max = lim_sup, out = enxame_fit_novo[:,0:-7]),enxame_estat),axis=1)   

        end = time.time()

        elapsed = end - start

        j.append(global_best[-1])

        obj.append(global_best[-7])

        pen_v.append(global_best[-6])

        pen_gq.append(global_best[-5])

        pen_tap.append(global_best[-4])

        pen_bsh.append(global_best[-3])
        
        pen_slack.append(global_best[-2])
        
        tempo.append(elapsed)
      

        # if relatorio == True:
            
        #     print(' ')

        #     print('Melhor Global da Iteração:',i)

        #     print('Objetivo:', global_best[-7])

        #     print('Penalização de Tensão:', global_best[-6])

        #     print('Penalização de Geração de Reativo:', global_best[-5])
            
        #     print('Penalização do Gerador Slack (P):', global_best[-2])

        #     print('Fitness:', global_best[-1])

            
        #     print('Tempo: ', elapsed)

        #     print(' ')

        #     print('_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ ')
            
            
    global_best[n_vgen:n_vgen+n_tap] = discreto_tap(global_best,n_tap,n_vgen,n_bshunt,sep)
        
    global_best[n_vgen+n_tap:n_vgen+n_tap+n_bshunt] = discreto_bshunt(global_best,n_tap,n_vgen,n_bshunt,sep)
        
    global_best[n_vgen+n_tap+n_bshunt:n_vgen+n_tap+n_bshunt+n_gens] = poz(n_vgen, n_tap, n_bshunt, n_gens, global_best, sep)[n_vgen+n_tap+n_bshunt:n_vgen+n_tap+n_bshunt+n_gens] 

    print(' ')

    print('Melhor Global da Iteração:',i)

    print('Objetivo:', global_best[-7])

    print('Penalização de Tensão:', global_best[-6])

    print('Penalização de Geração de Reativo:', global_best[-5])
            
    print('Penalização do Gerador Slack (P):', global_best[-2])

    print('Fitness:', global_best[-1])

            
    print('Tempo: ', elapsed)

    print(' ')

    print('_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ ')
                               
            
    return j,obj,pen_v,pen_gq,pen_tap,pen_bsh, pen_slack, global_best, tempo
