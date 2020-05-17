# -*- coding: utf-8 -*-

#%% Classe Das Partículas
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from scipy.integrate import odeint

class Particula:
    
    """
    
        Recebe no mínimo dois parâmetros para retornar todos os dados os dados da partícula descrita.
        Procura-se deixar o programa extremamente customizável para atender a magnetohidrodinâmica.
        >>= NECESSÀRIO inicializar Variáveis da partícula
        
        Parâmetros
        ----------
        vB: List EXCLUSIVAMENTE tridimensional, exemplo => (0,0,0)
            Campo de indução magnética do sistema.
        vE: List EXCLUSIVAMENTE tridimensional, exmplo => (0,0,0)
            Campo elétrico do sistema.
        n:  int (opcional),     Valor padrão => 100
            Número de subdivisões do Tempo T.
        T:  list (opcional) EXCLUSIVAMENTE bidimensional, Valor padrão => [0,1]
            Tempo em segundos do primeiro intervalo e do segundo.
        q:  float (opcional),   Valor padrão => 1.6e-19
            Valor da carga da partícula.
        m:  float (opcional),   Valor padrão => 9.11e-22
            Valor da massa da partícula.
        vi: list (opcional),   Valor padrão => [0.0, 0.0, 0.0]
            Valor da velocidade inicial da partícula.
        coordenadas: list (opcional),   Valor padrão => [0.0, 0.0, 0.0]
            Valor da coordenada inicial da partícula.
        
        Objetos
        -------
        Os objetos podem ser acessados para extrair ou modificar os dados.
        
        B:  np.Array tridimensional
            Retorna Vetor Campo de indução magnética. 
        E:  np.Array tridimensional
            Retorna Vetor Campo elétrico.
        n:  int
            Retorna Número de componentes do vos Vetores t, X, Y e Z
        t:  np.Array
            Retorna um Vetor das divisões do intervalo de tempo selecionado.
        q:  float
            Retorna a carga da partícula.
        m:  float
            Retorna a massa da partícula.
        w:  float
            Retorna a frequência de oscilação da partícula.
        r:  float
            Retorna o raio de Drift da partícula.
        X:  np.Array
            Retorna um Vetor das Divisões de espaços na coordena X do eixo. \n
            >>=[VARIÁVEL DEVE SER INICIALIZADA], Objeto => Posicoes3D
        Y:  np.Array
            Retorna um Vetor das Divisões de espaços na coordena Y do eixo. \n
            >>=[VARIÁVEL DEVE SER INICIALIZADA], Objeto => Posicoes3D
        Z:  np.Array
            Retorna um Vetor das Divisões de espaços na coordena Z do eixo. \n
            >>=[VARIÁVEL DEVE SER INICIALIZADA], Objeto => Posicoes3D
        vi: np.Array 
            Retorna o vetor velocidade inicial da partícula.
        mE: float
            Retorna o módulo do campo elétrico.
        mB: float
            Retorna o módulo do campo magnético de indução.
        vD: float
            Retorna a velocidade de Drift da partícula.
            
        
    """
    
    def __init__(self,vB:list,vE:list,n=1000,T=[0,1],q=1.6e-19,m=9.11e-22,vi=[0.0, 0.0, 0.0], coordenadas = [0.0 ,0.0 ,0.0]):
        
        assert type (vB) == list,   "VariavelERROR: O componente vB deve ser uma lista!"
        assert type (vE) == list,   "VariavelERROR: O componente vE deve ser uma lista!"
        assert len  (vB) == 3,      "VariavelERROR: O componente vB deve ser tridimensional!"
        assert len  (vE) == 3,      "VariavelERROR: O componente vE deve ser tridimensional!"
        assert len   (T) == 2,      "VariavelERROR: O componente T deve conter apenas 2 espaços!"
        
        self.B  = np.array(vB)
        self.E  = np.array(vE)
        self.n  = n
        self._T = T
        self.t  = np.linspace(self._T[0],self._T[1],self.n)  #Vetor com os tempos dos intervalos solicitados
        self.q  = q
        self.m  = m
        self.vi = np.array(vi)
        self.mE = np.linalg.norm(self.E)                    # Módulo de E
        self.mB = np.linalg.norm(self.B)                    # Módulo de B
        self.vD = self.mE/self.mB                           #Velocidade de Drift
        self.w  = self.q*self.mB/self.m                     # Frência de oscilação
        self.r  = self.m*self.vD/(self.q*self.mB)           # Raio de Drift
        self.X  = np.zeros(self.n)
        self.Y  = np.zeros(self.n)
        self.Z  = np.zeros(self.n)
        self.coordenadas = coordenadas
        
        
    def calcPosicao(self, export = False, nome = "Posicao3D"):
        """
        Calcula as posições da partícula, permite o retorno e exportação das triggers de dados
        
        Parâmetros
        ----------
        
        export: boll (Opcional)
                    Permite Salvar um arquivo .mhd (Magneto-HidroDinâmica). \n
                    >>= este arquivo contem os dados do rastreameno da partícula analisada detalhadamente
            nome:   str (Opcional)
                    Atribui um nome ao arquivo exportado pelo programa.
        """
        
        E  =     self.E
        B  =     self.B
        m  =     self.m
        q  =     self.q
        c  =     self.coordenadas
        v0 =     self.vi
        t  =     self.t
        c  =     self.coordenadas
        p0 =     [c[0],v0[0],c[1],v0[1],c[2],v0[2]]
        
        def edo(p,t):
            x,vx,y,vy,z,vz = p  
            dvx = q/m*(E[0]+ vy*B[2] - vz*B[1])
            dvy = q/m*(E[1]+ vz*B[0] - vx*B[2])
            dvz = q/m*(E[2]+ vx*B[1] - vy*B[0])
            s = [vx,dvx,vy,dvy,vz,dvz]
            return s
        
        ssol = odeint(edo, p0, t)
        
        u = ssol[:,1]
        v = ssol[:,3]
        w = ssol[:,5]
        
        self.X = ssol[:,0]
        self.Y = ssol[:,2]
        self.Z = ssol[:,4]
        self.V = np.array([ u, v, w ])
         
         #Exportar arquivo de dados:
        if export == True:
            hoje = datetime.now()
            arquivo = open(nome + '.mhd','a')
            arquivo.write('%%C  ' + hoje.strftime('%d/%m/%Y %H:%M') + '\n')
            i=0
            for i in range(self.n):
                arquivo.write(str(self.X[i]) + '  ||  ' + str(self.Y[i]) + '  ||  ' + str(self.Z[i]) + '\n')
            arquivo.close()
            
        return self.X, self.Y, self.Z
    
    def plot3D(self, Titulo='Plotagem3D', EixoX = 'X', EixoY = 'Y', EixoZ = 'Z',Salvar = False,Extencao = 'svg'):
        """
            Plota um gráfico 3D do deslocamento da partícula a cada instante
            
            Parâmetros
            ----------
            Titulo:     str (Opcional)
                        Atribui um título ao gráfico.
            EixoX:      str (Opcional)
                        Atribui um nome ao eixo X.
            EixoY:      str (Opcional)
                        Atribui um nome ao eixo Y.
            EixoZ:      str (Opcional)
                        Atribui um nome ao eixo Z.
            Salvar:     boll (Opcional)
                        Permite salvar um arquivo com o nome sendo igual ao titulo do gráfico.
            Extencao:   str (Opcional)
                        Atribui uma extenção ao arquivo salvo. \n
                        Por padrão é salvo como 'svg' mas aceita todos os outros formatos de imagem.
            
        """
        
        token = False
        for i in range(self.n):
            if not self.X[i] == 0 and self.Y[i] == 0 and self.Z[i] == 0:
                token = True
        
        #assert token == True, 'PlotagemERROR: Componentes X, Y e Z vazias, Execute Posicoes3D() antes de Plotar!!'
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(self.X, self.Y, self.Z, zdir='z')
        ax.set_xlabel(EixoX)
        ax.set_ylabel(EixoY)
        ax.set_zlabel(EixoZ)
        
        plt.title(Titulo)
        if Salvar == True:
            fig.savefig(Titulo + '.' + Extencao)
        plt.show()
        
    def plot2D(self, Titulo = 'Plotagem2D',EixoPadrao = 'xy', EixoX = 'X', EixoY = 'Y', Salvar = False, Extencao = 'svg'):
        """
            Plota um gráfico 3D do deslocamento da partícula a cada instante
            
            Parâmetros
            ----------
            Titulo:      str (Opcional)
                         Atribui um título ao gráfico.
            EixoPadrao:  str (Opcional)
                         Seleciona o eixo a ser utilisado. \n
                         >>= aceita os parâmetros:
                            ~~> 'xy' \n
                            ~~> 'xz' \n
                            ~~> 'yz' \n
                            
            EixoX:       str (Opcional)
                         Atribui um nome ao eixo X.
            EixoY:       str (Opcional)
                         Atribui um nome ao eixo Y.
            Salvar:      boll (Opcional)
                         Permite salvar um arquivo com o nome sendo igual ao titulo do gráfico.
            Extencao:    str (Opcional)
                         Atribui uma extenção ao arquivo salvo. \n
                         Por padrão é salvo como 'svg' mas aceita todos os outros formatos de imagem.
            
        """
        
        token = False
        for i in range(self.n):
            if not self.X[i] == 0 and self.Y[i] == 0 and self.Z[i] == 0:
                token = True
        
        #assert token == True, 'PlotagemERROR: Componentes X, Y e Z vazias, Execute Posicoes3D() antes de Plotar!!'
        
        EixoPadrao = EixoPadrao.lower()
        assert EixoPadrao == 'xy' or EixoPadrao == 'xz' or EixoPadrao == 'yz','PlotagemERROR: Variável EixoPadrao Preenchida incorretamente!!'
        
        fig, ax = plt.subplots()
        ax.set(xlabel = EixoX, ylabel = EixoY, title = Titulo)
        
        if EixoPadrao == 'xy' :
            ax.plot(self.X,self.Y)
        elif EixoPadrao == 'xz' :
            ax.plot(self.X,self.Z)
        elif EixoPadrao == 'yz' :
            ax.plot(self.Y,self.Z)
        
        if Salvar == True:
            fig.savefig(Titulo + '.' + Extencao)
        plt.show()
        
    def animate2D(self,Nome="Animacao2D",EixoPadrao = 'xy', Intervalo=5):
        """

        Parâmetros
        ----------
        Nome:           str (opcional)
                         Passa o nome do arquivo que será salvo.
        EixoPadrao:     str (opcional)
                         Seleciona o plano que será animado.\n
                            >>= aceita os parâmetros:
                                ~~> 'xy' \n
                                ~~> 'xz' \n
                                ~~> 'yz' \n
        Intervalo:      int (optional)
                         Seleciona o tempo em ms(milisegundos) entre os frames animados.

        """
        import matplotlib.animation as animate
        
        EixoPadrao = EixoPadrao.lower()
        assert EixoPadrao == 'xy' or EixoPadrao == 'xz' or EixoPadrao == 'yz','AnimateERROR: Variável EixoPadrao Preenchida incorretamente!!'
        
        fig, ax = plt.subplots()
        line, = ax.plot(0, 0)
        Max = [0,0]
        Min = [0,0]
        A   = []
        B   = []
        
        if EixoPadrao == 'xy' :
            #Mecanismo de Ordenamento dos Eixos
            MecX = self.X
            MecY = self.Y
            Max[0] = MecX[0]
            Max[1] = MecY[0]
            Min[0] = MecX[0]
            Min[1] = MecY[0]
            i=-1
            for i in range(self.n-1):
                if MecX[i] >= Max[0]:
                    Max[0] = MecX[i]
                if MecY[i] >= Max[1]:
                    Max[1] = MecY[i]
                if MecX[i] <= Min[0]:
                    Min[0] = MecX[i]
                if MecY[i] <= Min[1]:
                    Min[1] = MecY[i]
            Max[0] *= 1.1
            Max[1] *= 1.1
            Min[0] *= 1.1
            Min[1] *= 1.1
            ax.set_xlim(Min[0],Max[0])
            ax.set_ylim(Min[1],Max[1])
           #Fim do Mecanismo     
            def f(x):
                i=0
                while (x != self.X[i]):
                    i += 1
                B.append(self.Y[i-1])
                A.append(self.X[i-1])
                
                line.set_xdata(A)
                line.set_ydata(B)
                return line,
            animation = animate.FuncAnimation(fig, func = f,frames=self.X, interval = Intervalo )
        elif EixoPadrao == 'yz' :
            #Mecanismo de Ordenamento dos Eixos
            MecX = self.Y
            MecY = self.Z
            Max[0] = MecX[0]
            Max[1] = MecY[0]
            Min[0] = MecX[0]
            Min[1] = MecY[0]
            i=-1
            for i in range(self.n-1):
                if MecX[i] >= Max[0]:
                    Max[0] = MecX[i]
                if MecY[i] >= Max[1]:
                    Max[1] = MecY[i]
                if MecX[i] <= Min[0]:
                    Min[0] = MecX[i]
                if MecY[i] <= Min[1]:
                    Min[1] = MecY[i]
            Max[0] *= 1.1
            Max[1] *= 1.1
            Min[0] *= 1.1
            Min[1] *= 1.1
            ax.set_xlim(Min[0],Max[0])
            ax.set_ylim(Min[1],Max[1])
           #Fim do Mecanismo   
             
            def f(x):
                i=0
                while (x != self.Y[i]):
                    i += 1
                B.append(self.Z[i-1])
                A.append(self.Y[i-1])
                
                line.set_xdata(A)
                line.set_ydata(B)
                return line,
            animation = animate.FuncAnimation(fig, func = f,frames=self.Y, interval = Intervalo )
        elif EixoPadrao == 'xz' :
            #Mecanismo de Ordenamento dos Eixos
            MecX = self.X
            MecY = self.Z
            Max[0] = MecX[0]
            Max[1] = MecY[0]
            Min[0] = MecX[0]
            Min[1] = MecY[0]
            i=-1
            for i in range(self.n-1):
                if MecX[i] >= Max[0]:
                    Max[0] = MecX[i]
                if MecY[i] >= Max[1]:
                    Max[1] = MecY[i]
                if MecX[i] <= Min[0]:
                    Min[0] = MecX[i]
                if MecY[i] <= Min[1]:
                    Min[1] = MecY[i]
            Max[0] *= 1.1
            Max[1] *= 1.1
            Min[0] *= 1.1
            Min[1] *= 1.1
            ax.set_xlim(Min[0],Max[0])
            ax.set_ylim(Min[1],Max[1])
           #Fim do Mecanismo   
            
            def f(x):
                i=0
                while (x != self.X[i]):
                    i += 1
                B.append(self.Z[i-1])
                A.append(self.X[i-1])
                
                line.set_xdata(A)
                line.set_ydata(B)
                return line,
            animation = animate.FuncAnimation(fig, func = f,frames=self.X, interval = Intervalo )
            
        animation.save(filename = Nome + ".gif")
        plt.close()
        print("Animação Salvo com sucesso!")
        
        
        

