from resources.msc_rafael_pavan import *


def Comprehensive_Learning_PSO(sep, pt, rgap, zeta, psi, max_iter, n_particles, c1, v_amp, valor_inicial, step, wmin, wmax, relatorio=True,inicial=False):
        
    enxame_fit = cria_enxame(sep,n_particles)
            
    if len(sep.bus) == 14:
        
        n_vgen = 4+1
        n_tap = 3
        n_bshunt = 1
    
    if len(sep.bus) == 30:
        
        n_vgen = 5+1
        n_tap = 4
        n_bshunt = 2
        
    
    if len(sep.bus) == 118:
        
        n_vgen = 53+1
        n_tap = 9
        n_bshunt = 14
        
        
    if len(sep.bus) == 300:
        
        n_vgen = 68+1
        n_tap = 62
        n_bshunt = 29
    

    if inicial == True:
        
        enxame_fit[0,:]=valor_inicial   
        
    
    w_max=wmax
    
    w_min=wmin
    
    j = []
        
    tempo = []
        
    perdas = []
    
    pen_v = []
    
    pen_gq = []
    
    pen_tap = []
    
    pen_bsh = []

    refresh_rate = np.zeros((n_particles))
    
    v_lim_superior = np.repeat(sep.bus['max_vm_pu'][0], n_vgen)
    
    v_lim_inferior = np.repeat(sep.bus['min_vm_pu'][0], n_vgen)
    
    tap_pos, tap_neutral, tap_step_percent,valores_taps = coleta_dados_trafo(sep,relatorio=False)
    
    tap_max = np.repeat(valores_taps[-1], len(tap_pos))
    
    tap_min = np.repeat(valores_taps[0], len(tap_pos))
    
    bsh,b=coleta_dados_bshunt(sep)

    bsh_max=[]
    
    bsh_min=[]
    

    for bs in bsh:
    
        bsh_max.append([np.max(bs)])
    
        bsh_min.append([np.min(bs)])


    maximo = np.expand_dims(np.concatenate((v_lim_superior, tap_max, bsh_max), axis = None), 0)
    
    minimo = np.expand_dims(np.concatenate((v_lim_inferior, tap_min, bsh_min), axis = None), 0)
     
    
    lim_sup = np.tile(maximo, (n_particles,1))
    
    lim_inf = np.tile(minimo, (n_particles,1))
    
    delta = (lim_sup-lim_inf)[0,:]
    
    v_anterior = v_amp*cria_enxame(sep,n_particles)
    
    refresh_rate = np.zeros((n_particles))

    for i in range(0,max_iter):
       
            
        start = time.time()
        
        
        r1 = np.random.random_sample(size = (n_particles,1))
        
       
        enxame_fit_d = np.copy(enxame_fit)
    
        for linha in range(n_particles):
          
            enxame_fit_d[linha][n_vgen:n_vgen+n_tap] = discreto_tap(enxame_fit[linha].copy(),n_tap,n_vgen,n_bshunt,sep)
       
            enxame_fit_d[linha][n_vgen+n_tap:n_vgen+n_tap+n_bshunt] = discreto_bshunt(enxame_fit[linha].copy(),n_tap,n_vgen,n_bshunt,sep)
  

        
        enxame_fit[:,-6:] = (fluxo_de_pot(enxame_fit_d,sep))[:,-6:]
     
        enxame_fit[:,-6:] = (fitness(enxame_fit,zeta,psi,0,0))[:,-6:]
                       

        if i==0:
            
            best_particles = enxame_fit.copy()

            global_best = best_particles[np.argsort(best_particles[:, -1])][0,:].copy()
                
            pbestfi = np.copy(best_particles)
        

            
        for t in range(0,n_particles):
                
            if (enxame_fit[t,-1] < best_particles[t,-1]):
        
                best_particles[t,:] = enxame_fit[t,:].copy()

                refresh_rate[t] = 0

            else: refresh_rate[t] = refresh_rate[t] + 1

        
                            
        global_best = best_particles[np.copy(np.argsort(best_particles[:, -1]))][0,:].copy()

        
        for li in range(n_particles):
            
            if refresh_rate[li] >= rgap:
                
                pbestfi[li,:]=best_particles[li,:].copy()
                
                refresh_rate[li] = 0

                for co in range(len(global_best)-7):
                                          
                    if np.random.rand()<= pt:
                 
                        p1 = int(np.floor(n_particles*np.random.rand()))

                        p2 = int(np.floor(n_particles*np.random.rand()))

                        while(p1==p2):


                            p1 = int(np.floor(n_particles*np.random.rand()))

                            p2 = int(np.floor(n_particles*np.random.rand()))


                        part1 = best_particles[p1,:].copy()

                        part2 = best_particles[p2,:].copy()

                        if part1[-1] <= part2[-1]:

                            pbestfi[li,co] = part1[co].copy()

                        else:

                            pbestfi[li,co] = part2[co].copy()

                    
                flag = 0

                for p in range(len(pbestfi[li,:])-6):

                    if pbestfi[li,p] == best_particles[li,p]:

                        flag = flag+1

                if flag == len(pbestfi[li,:])-6:

                    l__ = int(np.floor(n_particles*np.random.rand()))

                    c__ = int(np.floor(len((global_best)-7)*np.random.rand()))
                           
                    part3 = best_particles[l__,:].copy()
                           
                    pbestfi[li,co] = part3[c__].copy()


        enxame_fit_anterior = enxame_fit.copy()
        
        w_novo = w_max-(w_max-w_min)*(i+1)/max_iter

        v_novo = np.multiply(w_novo,v_anterior.copy()) + c1*np.multiply(r1,(pbestfi-enxame_fit.copy())) 

        for l in range(v_novo.shape[0]):
                
                    for c in range(len(delta)):
                        
                        if v_novo[l,c] > delta[c]*step:
                            
                            v_novo[l,c] = delta[c].copy()*step
                            
                            
                        if v_novo[l,c] < -delta[c]*step:
                            
                            v_novo[l,c] = -delta[c].copy()*step
                            
                
                                                
        enxame_fit_novo = enxame_fit_anterior  + v_novo
        
        v_anterior = v_novo.copy()


        enxame_estat = enxame_fit_novo[:,-6:]

        enxame_fit = np.concatenate(( np.clip(enxame_fit_novo[:,0:-6], a_min = lim_inf, a_max = lim_sup, out = enxame_fit_novo[:,0:-6]),enxame_estat),axis=1)   
        
    
        end = time.time()

        elapsed = end - start

        j.append(global_best[-1])

        perdas.append(global_best[-6])

        pen_v.append(global_best[-5])

        pen_gq.append(global_best[-4])

        pen_tap.append(global_best[-3])

        pen_bsh.append(global_best[-2])
        
        tempo.append(elapsed)
      

        if relatorio == True:
            
            print(' ')

            print('Melhor Global da Iteração:',i)

            print('Perdas (pu):', global_best[-6])

            print('Penalização de Tensão:', global_best[-5])

            print('Penalização de Geração de Reativo:', global_best[-4])

            print('Fitness:', global_best[-1])
            
            print('Tempo: ', elapsed)

            print(' ')

            print('_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ ')
            

    global_best[n_vgen:n_vgen+n_tap] = discreto_tap(global_best,n_tap,n_vgen,n_bshunt,sep)
        
    global_best[n_vgen+n_tap:n_vgen+n_tap+n_bshunt] = discreto_bshunt(global_best,n_tap,n_vgen,n_bshunt,sep)
        

    return j, perdas, pen_v, pen_gq, pen_tap, pen_bsh, global_best, tempo
            