import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button,Slider
import easygui

#funcion para interpolar entre dos puntos
def lerp(t, x0, x1, y0, y1):
    return (1 - t)* x0 + t * x1, (1 - t)* y0 + t * y1

#Algoritmo de casteljau para generar curvas de bezier
def casteljau(x, y, t, all):
    n = len(x) - 1
    max_n = n
    nx = []
    ny = []
    nx.append(x)
    ny.append(y)

    while(n):
        ax =  []
        ay =  []
        for i in range(n):
            lp = lerp(t, nx[max_n - n][i], nx[max_n - n][i + 1], ny[max_n - n][i], ny[max_n - n][i + 1])
            ax.append(lp[0])
            ay.append(lp[1])
        n -= 1
        nx.append(ax)
        ny.append(ay)
    if (all == False):
        return nx[len(nx) - 1][0], ny[len(ny) - 1][0]
    return nx, ny


class InteractiveBezierMaker():
    epsilon = 5

    def __init__(self, fig, scatter, line):
        #Agregamos los botones y sus funciones
        ax1 = fig.add_axes([0.125, 0.9, 0.1, 0.075])
        saveButton = Button(ax1, 'Save')
        ax2 = fig.add_axes([0.24, 0.9, 0.1, 0.075])
        loadButton = Button(ax2, 'Load')
        ax3 = fig.add_axes([0.71, 0.9, 0.19, 0.075])
        sCasteljau = Button(ax3, 'Show Casteljau')
        saveButton.on_clicked(self.savePoints)
        loadButton.on_clicked(self.loadPoints)
        sCasteljau.on_clicked(self.showCasteljau)
        #self.show nos indica si debemos mostrar las lineas del casteljau o no
        #lo cambiamos a True o False al presionar el boton show casteljau
        self.show = False

        #Variables que nos sirven para la creacion de las curvas
        self.axis = ax #Eje de las rectas y las curvas
        self.scatter = scatter # Aqui gauaredamos la grafica de puntos y las contrucciones de bezier
        self.line = line
        self.lines = []
        self.colors = ['c', 'm', 'g', 'y', 'r'] #Colores que pueden tomar los diferentes niveles de bezier
        
        
        
        self.xs,  self.ys  = list(line.get_xdata()), list(line.get_ydata())  #Aqui guardamos las coordenadas x y y de los puntos de control
        self.px, self.py = [], [] #Guarda las coordenadas de la curva de bezier
        self._ind = None #indice del punto mas cercano al borrar puntos
        self.ax = scatter.axes #Los ejes donde est la grafica de puntos
        self.canvas = self.ax.figure.canvas # El canvas donde ploteamos las graficas
        #Al presionar el mouse se llama a la funcion indicada
        self.canvas.mpl_connect('button_press_event', self.button_press_callback)
        #Al soltar el mouse se llama a la funcion indicada
        self.canvas.mpl_connect('button_release_event', self.button_release_callback)
        #Al mover el mouse se llama a la funcion indicada
        self.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        #Al presionar una tecla se llama a la funcion indicada
        self.canvas.mpl_connect('key_press_event', self.on_key_press)
        
        #T es el valor que cambiaremos en el slider
        self.t = 0.0
        #definimos el slider
        axfreq = plt.axes([0.19, 0.0, 0.65, 0.03])
        self.freq_slider = Slider(
            ax=axfreq,
            label='casteljau(t)',
            valmin=0.0,
            valmax=1.0,
            valinit=self.t,
        )
        #al movernos en el slider actualizamos con la funcion update
        self.freq_slider.on_changed(self.update)
        #Ploteamos
        plt.show()    
    
    #funcion para guardar un conjunto de puntos
    def savePoints(self,event):
        path = easygui.fileopenbox(default = 'saved_points/')
        f = open(path, "w")
        n = len(self.xs)
        for i in range(n):
            f.write("%f\n" % self.xs[i])   
            f.write("%f\n" % self.ys[i])    
        f.close()
    #funcion para cargar puntos
    def loadPoints(self,event):
        #Abrimos la carpeta donde tenemos guardados puntos para cargar uno de ellos
        path = easygui.fileopenbox(default = 'saved_points/')
        #con el path obtenido abrimos el archivo
        f = open(path, "r")
        self.xs.clear()
        self.ys.clear()
        i = 1
        for line in f:
            if (i == 1): self.xs.append(float(line))
            if (i == -1): self.ys.append(float(line))
            i = -i
        mat = np.zeros((len(self.xs), 2))  
        mat[:,0] = self.xs
        mat[:,1] = self.ys
        scatter.set_offsets(mat)
        scatter.axes.figure.canvas.draw_idle()    
        self.updateGraph()

    #Muestra y borra el casteljau
    def showCasteljau(self,event): 
        if(self.show == False):
            self.show = True
            self.add_casteljau()
        elif(self.show == True):
            self.show = False
            self.erase_casteljau()
    #funcion para eliminar las lineas del casteljau    
    def erase_casteljau(self):
        first = False
        for i in self.axis.get_lines():
            if(first):
                i.remove()
            first = True
    #Funcion para crear la animacion del casteljau con un valor dado t
    def add_casteljau(self):
        bx, by = casteljau(self.xs, self.ys, self.t, True)
        for i in range(len(bx)):
            graph_color = self.colors[i % 5]
            self.ax.plot(bx[i], by[i], marker = 'o', color = graph_color)
    #Actualiza los niveles de la curva de bezier cuando se mueve el slider
    def update(self, t):
        self.t = t
        if(self.show):
            self.erase_casteljau()
            self.add_casteljau()
   
   #Actualiza la curva de bezier y su construccion cuando es necesario
    def updateGraph(self):
        self.px, self.py = [], []
        if(len(self.xs) > 1):
            for i in range(0, 1000):
                qtx, qty = casteljau(self.xs, self.ys, float(i/1000), False)
                self.px.append(round(qtx, 3))
                self.py.append(round(qty, 3))

            self.line.set_data(self.px, self.py)
                    
            if(self.show):
                self.erase_casteljau()
                self.add_casteljau()

        else:
            if(self.show):
                self.erase_casteljau()
            self.line.set_data([], [])

        self.line.figure.canvas.draw()

    #Funcion para ver si el punto de event ya esta en la lista de puntos(esto se ve si la distancia del punto
    # con respecto a uno de los de la lista es menor a eps)
    def get_ind_under_point(self, event):   
        if len(scatter.get_offsets()) == 0:
            return None
        xy = np.asarray(self.scatter.get_offsets())
        xyt = self.ax.transData.transform(xy)
        xt, yt = xyt[:, 0], xyt[:, 1]

        d = np.sqrt((xt - event.x)**2 + (yt - event.y)**2)
        ind = d.argmin()

        if d[ind] >= self.epsilon:
            ind = None

        return ind
    #Al presionar el mouse
    def button_press_callback(self, event):
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        #Obtenemos el indice del punto que se esta transladando    
        self._ind = self.get_ind_under_point(event)
    #Al soltar el mouse
    def button_release_callback(self, event):
        if event.button != 1:
            return
        self._ind = None
    
    #funcion para modificar la curva en tiempo real acorde a la translacion de uno de los puntos
    def motion_notify_callback(self, event):
        if self._ind is None:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        #Al indice que se esta moviendo modificamos sus coordenadas    
        x, y = event.xdata, event.ydata
        xy = np.asarray(self.scatter.get_offsets())
        xy[self._ind] = np.array([x, y])
        self.scatter.set_offsets(xy)
        self.canvas.draw_idle()
        
        self.xs[self._ind] = x
        self.ys[self._ind] = y

        self.updateGraph()

    #funcion que anade o elimina un punto y modifica la curva segun sea el cambio
    def on_key_press(self, event):
            if not event.inaxes:
                return
            #Si presionamos a entonces anadimos el punto en donde el mouse lo indica    
            if event.key == 'a' or event.key == 'A':
                self.xs.append(event.xdata)
                self.ys.append(event.ydata)
                new_point = [event.xdata, event.ydata]
                old_off = scatter.get_offsets()
                new_off = np.concatenate([old_off,np.array(new_point, ndmin=2)])
                scatter.set_offsets(new_off)
                scatter.axes.figure.canvas.draw_idle()
                self.updateGraph()

            #Si presionamos d, eliminamos el punto indicado por el mouse(de ser el caso)
            if event.key == 'd' or event.key == 'D':
                idx = self.get_ind_under_point(event)
                if idx is None:
                    return
    
                old_off = scatter.get_offsets()
                new_off = np.delete(old_off, idx, 0)
                scatter.set_offsets(new_off)
                scatter.axes.figure.canvas.draw_idle()

                self.xs.pop(idx)
                self.ys.pop(idx)
                

                self.updateGraph()

fig, ax = plt.subplots()
ax.set_title('Bezier Interactivo')
line, = ax.plot([], [])  # empty lin
scatter = ax.scatter([], [],  marker='o')


#Definimos los limites de la pantalla
plt.xlim(0, 10)
plt.ylim(0, 10)

#Inicializamos la clase
InteractiveBezierMaker(fig, scatter, line)