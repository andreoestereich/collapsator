#!/usr/bin/env python3

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure

import tkinter as tk
from tkinter import ttk
from tkinter.filedialog import askopenfilename
from tkinter import messagebox
from os.path import isfile
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import numpy as np

def reescale (exp, sizes, loaded_data):
     a = exp[0] # Y axis[2]
     b = exp[1] # X axis
     lamb = exp[2]
     reescaled = []
     #calculates the reescaled values
     for i in range(len(sizes)):
         reescaled.append([np.multiply(np.subtract(loaded_data[    i][0], lamb), pow(sizes[i], b)), np.multiply(loaded_data[i][1]    ,pow(sizes[i],a)), sizes[i]])
     return reescaled

def getLine(x,y):
    #makes a linear fit of the two neighbouring points in the reference curve
    slope = float(y[0] -  y[1])/(x[0] -  x[1])
    intercept = y[0] - slope* x[0]
    return np.poly1d([slope, intercept])

def distCurves(exp, sizes, loaded_data, xMin, xMax):
    reescaled = reescale(exp, sizes, loaded_data)
    yScaleFactor = np.power(max(sizes), exp[0])

    residual = 0.0

    #loops the reference curve over all curves
    for spreaded in range(len(loaded_data)):

        #loops over the remaining curves measuring the distance between them
        for idx in [x for x in range(len(loaded_data)) if x != spreaded]:
            #if xs[0].max > reescaled[spreaded][0].max():
            #    return 100000.0

            #loops over the points of the current curve
            for index in range(len(loaded_data[idx][0])):
                if loaded_data[idx][0][index] >= xMin and loaded_data[idx][0][index] <= xMax :

                    #finds the places where the reference curve is bigger than the current curve
                    wheres = np.where(reescaled[spreaded][0] > reescaled[idx][0][index])

                    #checks if there are points in the reference curve that are bigger than the current curve
                    if wheres[0].size != 0:
                        if wheres[0][0] != 0:

                            #makes a linear fit of the two neighbouring points in the reference curve
                            y = reescaled[spreaded][1][wheres[0][0] -1 : wheres[0][0]+1]
                            x = reescaled[spreaded][0][wheres[0][0] -1 : wheres[0][0]+1]
                            linha = getLine(x, y)
                            residual += abs(reescaled[idx][1][index] - linha(reescaled[idx][0][index]))/yScaleFactor

                        #if there are no points in the reference curve smaller than the current point it simply takes the distancec to the 1st point in the reference curve
                        else:
                            residual += abs(reescaled[idx][1][index] - reescaled[spreaded][1][0])/yScaleFactor

                    #in this case the reference curve is allways smaller than the current curve.
                    else:
                        residual += abs(reescaled[idx][1][index] - reescaled[spreaded][1][-1])/yScaleFactor
            #residual += np.sum(np.absolute(np.subtract(xs[1], fitted(xs[0]))))/len(xs[0])

    return residual


def validateFloat(action, index, value_if_allowed, prior_value, text, validation_type, trigger_type, widget_name):
    if text in '0123456789.-+':
        try:
            float(value_if_allowed)
            return True
        except ValueError:
            return False
    else:
        return False

def validateInt(action, index, value_if_allowed, prior_value, text, validation_type, trigger_type, widget_name):
    if text in '0123456789':
        try:
            int(value_if_allowed)
            return True
        except ValueError:
            return False
    else:
        return False

def validateFile(path):
    if not isfile(path):
        messagebox.showerror("Error", path + " is not a file.")
    return True

class inFrame:
    def __init__(self, frame):
        self.scrollCanvas = tk.Canvas(frame)
        self.frame = tk.Frame(self.scrollCanvas)
        self.scrollBar = tk.Scrollbar(frame, orient="vertical", command=self.scrollCanvas.yview)
        self.scrollBar.grid(column=1,row=0,sticky='nsew')
        self.scrollCanvas['yscrollcommand'] = self.scrollBar.set
        self.scrollBar.grid_forget()

        self.scrollCanvas.create_window((0,0),window=self.frame,anchor='nw')
        self.frame.bind("<Configure>", self.AuxscrollFunction)
        self.scrollCanvas.grid(column=0,row=0,sticky='nsew')

        addButton = ttk.Button(self.frame, text="Add file", command=self.addEntry)
        addButton.grid(column=1, row=100)
        #remButton = ttk.Button(self.frame, text="Remove file", command=self.delEntry)
        #remButton.grid(column=2, row=0)
        self.counter = 0
        self.addEntry()

    def addEntry(self):
        def cText(*args):
            entry.delete(0,'end')
            entry.insert(0,askopenfilename())

        def delEntry(*args):
            if self.counter > 1:
                entry.destroy()
                entryN.destroy()
                cBut.destroy()
                nEq.destroy()
                delBut.destroy()
                self.counter -= 1
        
        #this is important as to not violate the add button
        if self.counter < 99:
            self.counter+=1

            #to resize the scrollbar
            self.scrollBar.grid_forget()

            text = tk.StringVar()
            entry = ttk.Entry(self.frame, textvariable=text, width=30, validate = 'focusout', validatecommand = (self.frame.register(validateFile), '%s'))
            entry.grid(column=1, row=self.counter, columnspan=2)
            cBut = ttk.Button(self.frame, text="Open", command=cText)
            cBut.grid(column=0, row=self.counter, padx=(2,8))
            nEq = ttk.Label(self.frame, text=" N = ")
            nEq.grid(column=3, row=self.counter)

            vcmd = (self.frame.register(validateInt), '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')

            textN = tk.IntVar()
            textN.set(0)
            #textN.trace("w",callback=mustBeInt)
            entryN = ttk.Entry(self.frame, textvariable=textN, width=5, validate = 'key', validatecommand = vcmd)
            entryN.grid(column=4, row=self.counter)

            delBut = ttk.Button(self.frame, text="del", command=delEntry)
            delBut.grid(column=5, row=self.counter, padx=(8,2))

            self.scrollBar.grid(column=1,row=0,sticky='nsew')


    def fileLocations(self):
        entries = []
        for name in self.frame.winfo_children():
            if isinstance(name, ttk.Entry):
                entries.append(name.get())
        fileLocs = entries[0::2]
        sizes = entries[1::2]
        return (fileLocs, sizes)

    def enableChilds(self):
        for child in self.frame.winfo_children():
            if isinstance(child, ttk.Entry) or isinstance(child, ttk.Button):
                child.configure(state='normal')

    def disableChilds(self):
        for child in self.frame.winfo_children():
            if isinstance(child, ttk.Entry) or isinstance(child, ttk.Button):
                child.configure(state='disabled')
    
    def AuxscrollFunction(self,event):
        #You need to set a max size for frameTwo. Otherwise, it will grow as needed, and scrollbar do not act
        self.scrollCanvas.configure(scrollregion=self.scrollCanvas.bbox("all"),width=500,height=100)


class guessFrame:
    def __init__(self, frame):
        self.frame = frame

        ttk.Label(self.frame, text="a =").grid(row=0, column=0, sticky=(tk.E))
        ttk.Label(self.frame, text="b =").grid(row=1, column=0, sticky=(tk.E))
        ttk.Label(self.frame, text="x_c =").grid(row=2, column=0, sticky=(tk.E))
        ttk.Label(self.frame, text="x_min =").grid(row=3, column=0, sticky=(tk.E))
        ttk.Label(self.frame, text="x_max =").grid(row=4, column=0, sticky=(tk.E))

        vcmd = (self.frame.register(validateFloat), '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')

        self.aVar = tk.DoubleVar(value=0.)
        self.bVar = tk.DoubleVar(value=0.)
        self.xcVar = tk.DoubleVar(value=0.)
        self.xminVar = tk.DoubleVar()
        self.xmaxVar = tk.DoubleVar()

        aEntry = ttk.Entry(self.frame, textvariable=self.aVar, width=10, validate = 'key', validatecommand = vcmd, state='disabled')
        bEntry = ttk.Entry(self.frame, textvariable=self.bVar, width=10, validate = 'key', validatecommand = vcmd, state='disabled')
        xcEntry = ttk.Entry(self.frame, textvariable=self.xcVar, width=10, validate = 'key', validatecommand = vcmd, state='disabled')
        xminEntry = ttk.Entry(self.frame, textvariable=self.xminVar, width=10, validate = 'key', validatecommand = vcmd, state='disabled')
        xmaxEntry = ttk.Entry(self.frame, textvariable=self.xmaxVar, width=10, validate = 'key', validatecommand = vcmd, state='disabled')

        aEntry.grid(row=0,column=1)
        bEntry.grid(row=1,column=1)
        xcEntry.grid(row=2,column=1)
        xminEntry.grid(row=3,column=1)
        xmaxEntry.grid(row=4,column=1)

        #aVar.set(0.)
        #bVar.set(0.)
        #xcVar.set(0.)

        self.aFix = tk.BooleanVar()
        self.aFix.set(False)
        aFixButt = ttk.Checkbutton(self.frame, text='fixate', variable=self.aFix, state='disabled')
        #aFixButt.grid(row=0,column=2)

        self.bFix = tk.BooleanVar()
        self.bFix.set(False)
        bFixButt = ttk.Checkbutton(self.frame, text='fixate', variable=self.bFix, state='disabled')
        #bFixButt.grid(row=1,column=2)

        self.xcFix = tk.BooleanVar()
        self.xcFix.set(False)
        xcFixButt = ttk.Checkbutton(self.frame, text='fixate', variable=self.xcFix, state='disabled')
        #xcFixButt.grid(row=2,column=2)

    def enableChilds(self):
        for child in self.frame.winfo_children():
            if isinstance(child, ttk.Entry) or isinstance(child, ttk.Checkbutton):
                child.configure(state='normal')

    def disableChilds(self):
        for child in self.frame.winfo_children():
            if isinstance(child, ttk.Entry) or isinstance(child, ttk.Checkbutton):
                child.configure(state='disabled')

    def setDefault(self, guess, xmin, xmax):
        self.aVar.set(guess[0])
        self.bVar.set(guess[1])
        self.xcVar.set(guess[2])
        self.xminVar.set(xmin)
        self.xmaxVar.set(xmax)
    
    def returnVals(self):
        return [self.aVar.get(),self.bVar.get(),self.xcVar.get()]

class mainWin:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("Colapsator")

        #whole container
        self.content = ttk.Frame(self.root)

        #input files container
        filesFrame = ttk.LabelFrame(self.content, text="Input files", borderwidth=2)
        self.inFrame = inFrame(filesFrame)

        #plot preview container
        self.plotFrame = tk.Frame(self.content, borderwidth=5, relief="sunken")
        tk.Frame(self.plotFrame, width=300, height=300, bg='white').pack()

        self.distLabel = tk.Label(self.content, text="S = ????")

        #expoent and results frame
        guesFrame = ttk.Frame(self.content)
        self.resVal = guessFrame(guesFrame)

        self.Ninterac = tk.IntVar()
        self.Ninterac.set(100)
        self.interacEntry = tk.Entry(self.content, textvariable=self.Ninterac, validate = 'key', validatecommand = (self.content.register(validateInt), '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W'), width=8, state='disabled')
        self.approxButt = tk.Button(self.content, text="Approximate", command=self.approximate, state='disabled')
        
        self.statusLabel = tk.Label(self.content, text='')

        #adding things to the interface
        self.content.pack()
        filesFrame.grid(column=0, row=0, columnspan=3, padx=8, pady=8)
        self.plotFrame.grid(column=0, row=1, rowspan=4, padx=8, pady=4, sticky='nesw')

        self.loadedFlag = False
        self.loadButt = tk.Button(self.content, text="Load data", command=self.loadData)
        self.loadButt.grid(column=1, row=1, padx=8, pady=8)

        self.replotButt = tk.Button(self.content, text="Replot", state='disabled', command=self.graphReplot)
        self.replotButt.grid(column=2, row=1, padx=8, pady=8)

        guesFrame.grid(column=1, row=2, columnspan=2, padx=8, pady=4)
        self.interacEntry.grid(column=1, row=3, sticky=(tk.E, tk.S))
        tk.Label(self.content, text=" interactions").grid(column=2, row=3, sticky=(tk.W, tk.S))
        self.approxButt.grid(column=1, columnspan=2, row=4, sticky=(tk.N), pady=4)
        self.distLabel.grid(column=0, row=5, padx=4, pady=4)
        self.statusLabel.grid(column=1, columnspan=2, row=5, padx=4, pady=4)

        self.root.mainloop()

    def loadData(self):
        def loadFile(path):
            if isfile(path):
                try:
                    return np.loadtxt(path, usecols=[0,1])
                except ValueError:
                    return messagebox.askyesno("Error", path + " could not be loaded. Do you want to continue without it?")
            else:
                return messagebox.askyesno("Error", path + " is not a file. Do you want to continue without it?")

        if not self.loadedFlag:
            fileList, preSizes = self.inFrame.fileLocations()
            dataSet = []
            sizes = []

            for i in range(len(fileList)):
                res = loadFile(fileList[i])

                if isinstance(res, np.ndarray):
                    size= int(preSizes[i].replace('','0'))
                    if size== 0:
                        if not messagebox.askyesno("Error", "0 is not a valid system size. Do you want to continue without this file?"):
                            del dataSet[:]
                            del sizes[:]
                            return
                    else:
                        dataSet.append(res.transpose())
                        sizes.append(size)
                
                elif not res:
                    del dataSet[:]
                    del sizes[:]
                    return

            if len(dataSet) > 1:
                #print(dataSet[0][0])
                xmin = dataSet[0][0].min()
                xmax = dataSet[0][0].max()
                for xs in dataSet:
                    lmin = xs[0].min()
                    lmax = xs[0].max()
                    if lmin < xmin: xmin = lmin
                    if lmax > xmax: xmax = lmax

                self.resVal.setDefault([0.0,0.0,0.0],xmin,xmax)
                self.loadedFlag = True
                self.loadButt.configure(text="Change files")
                self.inFrame.disableChilds()
                self.resVal.enableChilds()
                self.replotButt.configure(state='normal')
                self.interacEntry.configure(state='normal')
                self.approxButt.configure(state='normal')

                self.plotData = dataSet
                self.sysSizes = sizes
                self.graphPlot()

            else:
                messagebox.showerror("Error", "You need at least 2 valid data sets.")

        else:
            self.loadedFlag = False
            self.loadButt.configure(text="Load data")
            self.resVal.disableChilds()
            self.inFrame.enableChilds()
            self.replotButt.configure(state='disabled')
            self.interacEntry.configure(state='disabled')
            self.approxButt.configure(state='disabled')

    def graphPlot(self):

        #update distance label
        self.updateDist()

        for child in self.plotFrame.winfo_children():
            child.destroy()

        f = Figure(figsize=(4,4), dpi=100)
        self.axis = f.add_subplot(111)

        for graph in self.plotData:
            self.axis.plot(graph[0],graph[1], 'o')

        self.canvas = FigureCanvasTkAgg(f, self.plotFrame)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(self.canvas, self.plotFrame)
        toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def graphReplot(self):
        #update info labels
        self.updateDist()
        self.statusLabel.configure(text='')

        reescaled = reescale(self.resVal.returnVals(),self.sysSizes, self.plotData)
        xmin = reescaled[0][0].min()
        xmax = reescaled[0][0].max()
        ymin = reescaled[0][1].min()
        ymax = reescaled[0][1].max()
        for xs in reescaled:
            lmin = xs[0].min()
            lmax = xs[0].max()
            if lmin < xmin: xmin = lmin
            if lmax > xmax: xmax = lmax
            lmin = xs[1].min()
            lmax = xs[1].max()
            if lmin < ymin: ymin = lmin
            if lmax > ymax: ymax = lmax
        
        self.axis.clear()
        self.axis.set_xlim(xmin, xmax)
        self.axis.set_ylim(ymin, ymax)

        for graph in reescaled:
            self.axis.plot(graph[0],graph[1], 'o')

        self.canvas.draw()

    def updateDist(self):
        dist = distCurves(self.resVal.returnVals(),self.sysSizes, self.plotData, self.resVal.xminVar.get(), self.resVal.xmaxVar.get())
        self.distLabel.configure(text="S = %g"%(dist))

    def approximate(self):
        guess = self.resVal.returnVals()
        xmin = self.resVal.xminVar.get()
        xmax = self.resVal.xmaxVar.get()
        
        res = minimize(distCurves, guess , args=(self.sysSizes, self.plotData, xmin, xmax), method='Nelder-Mead', options={'maxiter': self.Ninterac.get()})
        
        self.resVal.setDefault(res['x'],xmin,xmax)
        #self.distLabel.configure(text="S = %g"%(res['fun']))
        self.graphReplot()
        self.statusLabel.configure(text=res['message'])

mainWin()

