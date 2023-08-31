from tkinter import *
import tkinter.messagebox as tmsg
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import matplotlib.figure
import numpy as np



###################################################################################################
#the mathematical engine

x_min, x_max = -2, 2
num_points = 1000
x = np.linspace(x_min, x_max, num_points)
dx = (x_max-x_min)/num_points
zero = 0*x
length = x_max- x_min

t = np.linspace(0,1.5,500)

def non_phys_exp(x,t):
    return np.exp(t-x)

def non_phys_exp_2(x,t):
    return np.exp(x-t)

def phys_exp(x,t):
    a = ((x - t) - x_min) % length + x_min 
    return np.exp(-np.abs(a))

def L2_norm(psi, dx):
    norm = np.sqrt(np.trapz(np.abs(psi) ** 2, dx=dx))
    return norm

################################################################################################
class GUI(Tk):
    def __init__(self):
        super().__init__()
        #this fixes the size
        self.geometry("1600x700")
        self.maxsize(1600,700)
        self.minsize(1600,700)
        #title of our gui
        self.title("Boundary Conditions and Dynamics")

        self.choice = StringVar()
        self.choice.set("Radio")
        self.choice.set("e^x")

#########################################################################################################   
    def compute(self,*args):
        self.ax.cla()

        if self.choice.get()=="e^x":
            for i in range(0,len(t)):
                if i%10==0:
                    self.ax.plot(x, non_phys_exp_2(x,t[i]), c='blue')        
                    self.ax.set_title('Time Evolution of Chosen Initial Data')
                    self.ax.axvline(x=-2, c='white', linestyle='dashed')
                    self.ax.axvline(x=+2, c='white', linestyle='dashed')
                    self.ax.set_ylabel(r'$\psi(x,t)$')
                    self.ax.set_xlabel('Domain')    
                    self.ax.set_ylim([-0.1,1.1])             
                    self.canvas.draw()
                    self.canvas.start_event_loop(0.1)
                    self.ax.cla()
        
        elif self.choice.get()=="e^-x":
            for i in range(0,len(t)):
                if i%10==0:
                    self.ax.plot(x, non_phys_exp(x,t[i]), c='blue')        
                    self.ax.set_title('Time Evolution of Chosen Initial Data')
                    self.ax.axvline(x=-2, c='white', linestyle='dashed')
                    self.ax.axvline(x=+2, c='white', linestyle='dashed') 
                    self.ax.set_ylabel(r'$\psi(x,t)$') 
                    self.ax.set_xlabel('Domain') 
                    self.ax.set_ylim([-0.1,1.1])                 
                    self.canvas.draw()
                    self.canvas.start_event_loop(0.1)
                    self.ax.cla()

        elif self.choice.get()=="e^-|x|":
            for i in range(0,len(t)):
                if i%10==0:
                    self.ax.plot(x, phys_exp(x,t[i]), c='blue')        
                    self.ax.set_title('Time Evolution of Chosen Initial Data') 
                    self.ax.axvline(x=-2, c='white', linestyle='dashed')
                    self.ax.axvline(x=+2, c='white', linestyle='dashed')
                    self.ax.set_ylabel(r'$\psi(x,t)$')  
                    self.ax.set_xlabel('Domain')
                    self.ax.set_ylim([-0.1,1.1])                
                    self.canvas.draw()
                    self.canvas.start_event_loop(0.1)
                    self.ax.cla()
        
#########################################################################################################           
    def show_help(self):
        tmsg.showinfo("Tips", "This applet shows the time evolution of a wave function that is initially an exponential. The Hamiltonian is now just equal to momentum, so we observe translation of the wave function in time.")
#############################################################################################################    
    def make_title(self):
        f1 = Frame(self, bg='red', borderwidth=10, relief=GROOVE)
        f1.pack(side=TOP, fill=X, pady=10)
        Label(f1, bg = "red", fg = "black", font = ("calibiri", 20, "bold"),  text = "Boundary Conditions and Dynamics").pack()
##########################################################################################################################    
    def make_selection_frame(self):
        f2 = Frame(self, bg = "lightblue", borderwidth=10 , relief = RAISED)
        f2.pack(fill=X)
        Label(f2, bg = "lightblue", fg = "black", font = ("calibiri", 15, "bold"), text = "Initial Data Selection").pack()

        e_x = Radiobutton(f2, bg = "lightblue", fg='black',text="e^x",  variable=self.choice, value="e^x")
        e_x.pack(side=LEFT, padx=10)

        e_mx = Radiobutton(f2, bg = "lightblue", fg='black',text="e^-x",  variable=self.choice, value="e^-x")
        e_mx.pack(side=LEFT,padx=10)

        e_mmx = Radiobutton(f2, bg = "lightblue", fg='black',text="e^-|x|",  variable=self.choice, value="e^-|x|")
        e_mmx.pack(side=LEFT,padx=10)

        e_x.config(command=self.compute)
        e_mx.config(command=self.compute)
        e_mmx.config(command = self.compute)
        

        
##############################################################################################################################################
    def make_plot(self):
        f4 = Frame(self, bg = "lightblue", borderwidth = 10, relief = RAISED)
        f4.pack(fill=X)
        self.fig = matplotlib.figure.Figure(edgecolor='black', facecolor='lightblue',linewidth=7)
        self.fig.set_size_inches(13,4.5)
        self.ax = self.fig.add_subplot(1,1,1)

        self.ax.set_facecolor('#FF8D33')
        self.canvas = FigureCanvasTkAgg(self.fig, master=f4)  
        self.canvas.get_tk_widget().pack(pady=5)

        # Create Toolbar
        toolbar = NavigationToolbar2Tk(self.canvas, f4, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(fill=X)


        help = Button(f4, fg='black', bg='red', text="Help", font = ("calibiri", 12, "bold"), command = self.show_help)
        help.pack()
#############################################################################################################################################

#The main body of the program, which just calls all the functions from before

#instance of out GUI class
bcd = GUI() #stands for boundary conditions and dynamics

#methods that give life to the GUI
bcd.make_title()
bcd.make_selection_frame()
bcd.make_plot()
bcd.compute()

#the infinite GUI loop
bcd.mainloop()