from tkinter import *
import tkinter.messagebox as tmsg
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import matplotlib.figure
import numpy as np



###################################################################################################
#the mathematical engine

def create_tridiagonal_matrix(d, case):
    
    if case == 'real':
        # Generate random diagonal elements from a Gaussian distribution
        diagonal = np.random.normal(loc=0, scale=1, size=d)

        # Ensure the average is zero
        mean_diagonal = np.mean(diagonal)
        diagonal -= mean_diagonal

        # Ensure the variance is one
        variance_diagonal = np.var(diagonal)
        diagonal /= np.sqrt(variance_diagonal)
        
    elif case =='zero':
        # Diagonal elements are just zeros of length equal to d
        diagonal = np.zeros(d)
        
    
    # Create the tridiagonal matrix
    matrix = np.diag(diagonal)
    matrix_upper = np.ones(d - 1)
    matrix_lower = np.ones(d - 1)
    np.fill_diagonal(matrix[1:], matrix_upper)
    np.fill_diagonal(matrix[:, 1:], matrix_lower)

    #periodic boundary conditions
    matrix[0,d-1] = 1.0
    matrix[d-1, 0] = 1.0
    
    return matrix



################################################################################################
class GUI(Tk):
    def __init__(self):
        super().__init__()
        #this fixes the size
        self.geometry("1600x750")
        self.maxsize(1600,750)
        self.minsize(1600,750)
        #title of our gui
        self.title("Discrete Hamiltonian Spectrum")

        self.choice = StringVar()
        self.choice.set("Radio")
        self.choice.set("zero")
#####################################################################################################
    def change_plot(self, *args):
        self.ax.cla()
        self.ax2.cla()

        #eigen value distribution on the left
        self.ax.plot(self.eigenvalues, c='blue')
        self.ax.axvline(x=self.ev.get(), c='black', linestyle='dashed')
        self.ax.axhline(y=self.eigenvalues[self.ev.get()-1], c='black', linestyle='dashed')
        self.ax.set_title('Eigen Value Distribution')
        self.ax.set_xlabel('List Index')
        self.ax.set_ylabel('Eigen Values')
        
        #eigen vector plot on the right
        self.evec = self.eigenvectors[:,self.ev.get()-1]
        self.ax2.plot((self.evec)**2.0, c='blue')
        self.ax2.set_title('Corresponding Eigen Vector')
        self.ax2.set_xlabel('List Index')
        self.ax2.set_ylabel(r'$\left|\psi\right|^2$')
        self.ax2.set_ylim([np.min(self.evec**2.0)-0.001, np.max(self.evec**2.0)+0.001])

        self.canvas.draw()

 #########################################################################################################   
    def compute(self,*args):
        self.ev.config(from_=1, to=self.d.get())
        self.ax.cla()
        self.ax2.cla()

        # Create the tridiagonal matrix
        hamiltonian = create_tridiagonal_matrix(self.d.get(), self.choice.get())

        # Compute eigenvalues and eigenvectors
        self.eigenvalues, self.eigenvectors = np.linalg.eigh(hamiltonian)

        #eigen value distribution on the left
        self.ax.plot(self.eigenvalues, c='blue')
        self.ax.axvline(x=self.ev.get(), c='black', linestyle='dashed')
        self.ax.axhline(y=self.eigenvalues[self.ev.get()-1], c='black', linestyle='dashed')
        self.ax.set_title('Eigen Value Distribution')
        self.ax.set_xlabel('List Index')
        self.ax.set_ylabel('Eigen Values')
       

        #eigen vector plot on the right
        self.evec = self.eigenvectors[:,self.ev.get()-1]
        self.ax2.plot((self.evec)**2.0, c='blue')
        self.ax2.set_title('Corresponding Eigen Vector')
        self.ax2.set_xlabel('List Index')
        self.ax2.set_ylabel(r'$\left|\psi\right|^2$')
        self.ax2.set_ylim([np.min(self.evec**2.0)-0.001, np.max(self.evec**2.0)+0.001])

        self.canvas.draw()

#########################################################################################################           
    def show_help(self):
        tmsg.showinfo("Tips", "This applet shows the eigen values and eigen vectors of a random discrete Hamiltonian. Change the potential and dimension to see how the spectrum looks like.")
#############################################################################################################    
    def make_title(self):
        f1 = Frame(self, bg='red', borderwidth=10, relief=GROOVE)
        f1.pack(side=TOP, fill=X, pady=10)
        Label(f1, bg = "red", fg = "black", font = ("calibiri", 20, "bold"),  text = "Discrete Hamiltonian Spectrum").pack()
##########################################################################################################################    
    def make_selection_frame(self):
        f2 = Frame(self, bg = "lightblue", borderwidth=10 , relief = RAISED)
        f2.pack(fill=X)
        Label(f2, bg = "lightblue", fg = "black", font = ("calibiri", 15, "bold"), text = "Hamiltonian Selection").pack()

        v0 = Radiobutton(f2, bg = "lightblue", fg='black',text="V=0",  variable=self.choice, value="zero")
        v0.pack(side=LEFT, padx=10)

        vn = Radiobutton(f2, bg = "lightblue", fg='black',text="Random V",  variable=self.choice, value="real")
        vn.pack(side=LEFT,padx=10)

       

        Label(f2, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="Dimension").pack(side = LEFT, padx=10)
        self.d = Scale(f2, bg ='lightblue', fg='black',from_=20, to=200, orient=HORIZONTAL)
        self.d.set(100)
        self.d.pack(side = LEFT, padx=10)

        Label(f2, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="Eigen Value (Choose By Index)").pack(side = LEFT, padx=10)
        self.ev = Scale(f2, bg ='lightblue', fg='black',from_=1, to=self.d.get(), orient=HORIZONTAL)
        self.ev.set(1)
        self.ev.pack(side = LEFT, padx=10)

        v0.config(command=self.compute)
        vn.config(command=self.compute)
        self.d.config(command = self.compute)
        self.ev.config(command=self.change_plot)

        
##############################################################################################################################################
    def make_plot(self):
        f4 = Frame(self, bg = "lightblue", borderwidth = 10, relief = RAISED)
        f4.pack(fill=X)
        self.fig = matplotlib.figure.Figure(edgecolor='black', facecolor='lightblue',linewidth=7)
        self.fig.set_size_inches(13,4.5)
        self.ax = self.fig.add_subplot(1,2,1)
        self.ax2 = self.fig.add_subplot(1,2,2)

        self.ax.set_facecolor('#FF8D33')
        self.ax2.set_facecolor('#FF8D33')
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
discrete_hamiltonian = GUI()

#methods that give life to the GUI
discrete_hamiltonian.make_title()
discrete_hamiltonian.make_selection_frame()
discrete_hamiltonian.make_plot()
discrete_hamiltonian.compute()

#the infinite GUI loop
discrete_hamiltonian.mainloop()
