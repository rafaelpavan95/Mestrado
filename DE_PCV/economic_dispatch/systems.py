import numpy as np


class system_3():

    def __init__(self):
                        
        self.pg_max = np.array([600, 200, 400])
        self.pg_min = np.array([100, 50, 100])
        self.demanda = 850
        self.a_k = np.array([0.001562, 0.00482,0.00194])
        self.b_k = np.array([7.92, 7.97, 7.85])
        self.c_k = np.array([561, 78, 310])
        self.e_k = np.array([300, 150, 200])
        self.f_k = np.array([0.0315, 0.063, 0.042])
            
            
    def costs(self, potencias):
        
        return np.sum(self.a_k*potencias**2 + self.b_k*potencias + self.c_k + np.abs(self.e_k*np.sin(self.f_k*(self.pg_min-potencias))))
    
    
class system_13():

    def __init__(self):
            
        # Referência: João Vitor Dias, Diego Silva, Pothya

        self.pg_max = np.array([680, 360, 360, 180, 180, 180, 180, 180, 180, 120, 120, 120, 120]);
        self.pg_min = np.array([0, 0, 0, 60, 60, 60, 60, 60, 60, 40, 40, 55, 55]);
        self.demanda = 2520;
        self.a_k = np.array([0.00028, 0.00056, 0.00056, 0.00324, 0.00324, 0.00324, 0.00324, 0.00324, 0.00324, 0.00284, 0.00284, 0.00284, 0.00284]);
        self.b_k = np.array([8.1, 8.1, 8.1, 7.74, 7.74, 7.74, 7.74, 7.74, 7.74, 8.6, 8.6, 8.6, 8.6]);
        self.c_k = np.array([550, 309, 307, 240, 240, 240, 240, 240, 240, 126, 126, 126, 126]);
        self.e_k = np.array([300, 200, 200, 150, 150, 150, 150, 150, 150, 100, 100, 100, 100]);
        self.f_k = np.array([0.035, 0.042, 0.042, 0.063, 0.063, 0.063, 0.063, 0.063, 0.063, 0.084, 0.084, 0.084, 0.084]);

            
    def costs(self, potencias):
        
        return np.sum(self.a_k*potencias**2 + self.b_k*potencias + self.c_k + np.abs(self.e_k*np.sin(self.f_k*(self.pg_min-potencias))))
    
    
class system_19():

    def __init__(self):
        
        # Referência: Diego Nunes Silva, João Vitor Dias
        self.pg_min = np.array([100, 120, 100, 8, 50, 150, 50, 100, 200, 15, 50, 25, 50, 0, 20, 15, 15, 50, 400]);  
        self.pg_max = np.array([300, 438, 250, 25, 63.75, 300, 63.75, 500, 600, 40, 150, 75, 63.75, 95, 220, 80, 80, 230, 500]);
        self.a_k = np.array([0.0097, 0.0055, 0.0055, 0.0025, 0, 0.008, 0, 0.0075, 0.0085, 0.009, 0.0045, 0.0025, 0, 0.0045, 0.0065, 0.0045, 0.0025, 0.0045, 0.008]);
        self.b_k = np.array([6.8, 4, 4, 0.85, 5.28, 3.5, 5.439, 6, 6, 5.2, 1.6, 0.85, 2.55, 1.6, 4.7, 1.4, 0.85, 1.6, 5.5]);
        self.c_k = np.array([119, 90, 45, 0, 0.891, 110, 21, 88, 55, 90, 65, 78, 49, 85, 80, 90, 10, 25, 90]);
        self.e_k = np.array([90, 79, 0, 0, 0, 0, 0, 50, 0, 0, 0, 58, 0, 0, 92, 0, 0, 0, 0]); 
        self.f_k = np.array([0.72, 0.05, 0, 0, 0, 0, 0, 0.52, 0, 0, 0, 0.02, 0, 0, 0.75, 0, 0, 0, 0]); 
        self.demanda = 2908
            
    def costs(self, potencias):
        
        return np.sum(self.a_k*potencias**2 + self.b_k*potencias + self.c_k + np.abs(self.e_k*np.sin(self.f_k*(self.pg_min-potencias))))
    
    
class system_40():

    def __init__(self):
        
        # Referência: https://sci-hub.se/https://ieeexplore.ieee.org/document/1179910

        self.pg_min = np.array([36, 36,60,80,47,68,110,135,135,130,94,94,125,125, 125,125,220,220,242,242,254,254,254,254,254,254,10,10,10,47,60 ,60, 60, 90, 90, 90, 25, 25, 25, 242]);
        self.pg_max = np.array([114,114,120,190,97,140,300,300,300,300,375,375,500, 500,500,500,500,500,550,550,550,550,550, 550,550,550,150,150,150,97,190 ,190, 190, 200, 200, 200, 110, 110, 110, 550]);
        self.c_k = np.array([94.705,94.705,309.54,369.03,148.89,222.33,287.71,391.98,455.76,722.82,635.2,654.69,913.4,1760.4,1728.3,1728.3,647.85,649.69,647.83,647.81,785.96,785.96,794.53,794.53,801.32,801.32,1055.1,1055.1,1055.1,148.89,222.92, 222.920, 222.920, 107.870, 116.580, 116.580, 307.450, 307.450, 307.450, 647.830]);
        self.b_k = np.array([6.73,6.73,7.07,8.18,5.35,8.05,8.03,6.99,6.6,12.9,12.9,12.8,12.5,8.84,9.15,9.15,7.97,7.95,7.97,7.97,6.63,6.63,6.66,6.66,7.1,7.1,3.33,3.33,3.33,5.35,6.43, 6.43, 6.43, 8.95, 8.62, 8.62, 5.88, 5.88, 5.88, 7.97]);
        self.a_k = np.array([0.0069,0.0069,0.02028,0.00942,0.0114,0.01142,0.00357,0.00492,0.00573,0.00605,0.00515,0.00569,0.00421,0.00752,0.00708,0.00708,0.00313,0.00313,0.00313,0.00313,0.00298,0.00298,0.00284,0.00284,0.00277,0.00277,0.52124,0.52124,0.52124,0.0114,0.0016, 0.00160, 0.00160, 0.00010, 0.00010, 0.00010, 0.01610, 0.01610, 0.01610, 0.00313]);
        self.e_k = np.array([100, 100,100,150,120,100,200,200,200,200,200,200,300,300,300,300, 300,300,300,300,300,300,300,300,300,300,120,120,120,120,150, 150, 150, 200, 200, 200, 80, 80, 80, 300]);                         
        self.f_k = np.array([0.084, 0.084, 0.084, 0.063, 0.077, 0.084, 0.042, 0.042, 0.042, 0.042, 0.042, 0.042, 0.035, 0.035, 0.035, 0.035, 0.035, 0.035, 0.035, 0.035, 0.035, 0.035, 0.035, 0.035, 0.035, 0.035, 0.077, 0.077, 0.077, 0.077, 0.063, 0.063, 0.063, 0.042, 0.042, 0.042, 0.098, 0.098, 0.098, 0.035]);
        self.demanda = 10500
        
    def costs(self, potencias):
        
        return np.sum(self.a_k*potencias**2 + self.b_k*potencias + self.c_k + np.abs(self.e_k*np.sin(self.f_k*(self.pg_min-potencias))))
    
    
