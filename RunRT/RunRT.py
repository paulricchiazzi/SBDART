import sys, os
if len(sys.argv) == 1:
    import matplotlib
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D
    from Tkinter import Tk, Button, Frame, Spinbox, Menu, IntVar, BooleanVar, StringVar, Toplevel, Text, Entry, Checkbutton, Label
    import tkFileDialog
    import Ephemeris
    import Spinner
    import inspect
    import glob
    import datetime
    import pickle
import numpy as np
import platform
from collections import OrderedDict
import timeit
import GenInput
import RtReader
import subprocess
import os

sbdartexe = '/Users/paulricchiazzi/Documents/SBDART/sbdart'

class RunRT:
    """base class for radiative transfer GUI"""

    def __init__(self, master):

        # properties

        self.sbdartexe = sbdartexe

        self.geninput = GenInput.GenInput()     # instance of GenInput class
        self.rtdocwindow = None                 # help popup window object
        self.rtdoctext = None                   # help text object
        self.runname = "runrt"                  # run name
        self.sbdartoutput = ""                  # sbdart output text buffer
        self.variant_to_plot = 0                # variant parameter included on plot
        self.groupmenuSpinners = []             # plot group spin boxes
        self.optionSpinners = []                # sza to hour option spin boxes
        self.plotwindow=[]                      # xmin,ymin,xmax,ymax
        self.colorbarzoom=[]                    # normalized zmin,zmax,dzmin,dzmax
        self.plottype = 'xy'
        self.yrange = [float('inf'), float('-inf')]
        self.FixRangeStep = 0
        self.abortit = False
        self.rectXY = (0.125,0.1,0.9,0.9)
        self.rectPolar = (0.01,0.01,0.9,0.99)
        self.leg = None
        self.diffbase = ''
        self.rightmouse = "<Button-3>" if platform.system() == "Darwin" else "<Button-2>"  # OSX reverses Button-2 and Button-3

        # end properties

        self.master = master
        master.title(self.runname)

        self.menubar = Menu(master)
        self.filemenu = Menu(self.menubar, tearoff=0)
        self.filemenu.add_command(label="Recover", command=lambda: self.LoadFile("RUNS/~RunRT.sbd"))
        self.filemenu.add_command(label="Open", command=self.LoadFile)
        self.filemenu.add_command(label="Write", command=self.WriteFile)
        self.filemenu.add_command(label="Save as", command = self.PickleSave)
        self.filemenu.add_command(label="Save", command = lambda sfile = self.runname : self.PickleSave(sfile))
        self.filemenu.add_separator()
        #self.filemenu.add_command(label="View Output", command=self.ViewOutput)
        self.filemenu.add_command(label="Save plot as PNG", command = lambda choice = 'Save': self.PlotData(choice))
        self.filemenu.add_command(label="View Plot Data", command = lambda choice = 'View': self.PlotData(choice))
        #self.filemenu.add_command(label="Copy Plot Data", command = lambda choice = 'Copy': self.PlotData(choice))
        self.menubar.add_cascade(label="File",menu=self.filemenu)

        self.optionmenu = Menu(self.menubar, tearoff=0)
        self.optionDiurnalPlot = IntVar(0)
        self.optionmenu.add_radiobutton(label="Nominal", variable=self.optionDiurnalPlot, value=0, command=self.SetupNominalPlots)
        self.optionmenu.add_radiobutton(label='Hourly flux', variable=self.optionDiurnalPlot, value=1, command=self.SetupHourly)
        self.optionmenu.add_radiobutton(label='Hourly flux (fine)', variable=self.optionDiurnalPlot, value=2, command=self.SetupHourlyPlace)
        self.optionmenu.add_radiobutton(label='Diurnal average vs lat', variable=self.optionDiurnalPlot, value=3, command=lambda parm='day' : self.SetupDailySpinner(parm))
        self.optionmenu.add_radiobutton(label='Diurnal average vs day', variable=self.optionDiurnalPlot, value=4, command=lambda parm='lat' : self.SetupDailySpinner(parm))
        self.optionmenu.add_separator()
        self.optionComparisonPlot = IntVar(0)
        self.optionmenu.add_radiobutton(label='No Comparison', variable = self.optionComparisonPlot, value=0, command = self.Plotit)
        self.optionmenu.add_radiobutton(label='Difference Plot', variable = self.optionComparisonPlot, value=1, command = self.Plotit)
        self.optionmenu.add_radiobutton(label='Ratio Plot', variable = self.optionComparisonPlot, value=2,command = self.Plotit)
        self.menubar.add_cascade(label="Options", menu=self.optionmenu)

        # select plot groups

        self.menugroups = os.listdir("RUNS")
        self.menugroup = IntVar(value=-1)
        if len(self.menugroups) > 0:
            self.groupmenu = Menu(self.menubar, tearoff=0)
            if len(self.menugroups) == 0:
                self.groupmenu.add_radiobutton(label="None", variable=self.menugroup, value=-1)
            for k in range(0, len(self.menugroups)):
                self.groupmenu.add_radiobutton(label=self.menugroups[k], variable=self.menugroup, value=k,
                                               command=self.SetCaption)
                # self.menubar.add_cascade(label="Runs", menu=self.groupmenu)


            # setup menu for loop selection

        self.groupmenu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Plot Lines", menu=self.groupmenu)

        # setup menu for flux parameters

        self.fmenu=Menu(self.menubar, tearoff=0)
        self.parser = RtReader.RtReader("IOUT=10", BooleanVar)
        self.SetRtMenu()
        self.menubar.add_cascade(label="RT Quants", menu=self.fmenu)

        # help menu

        self.helpmenu = Menu(self.menubar, tearoff=0)
        self.helpmenu.add_command(label="RunRT", command=self.HelpGUI)
        self.helpmenu.add_command(label="RTdoc", command=self.HelpRtDoc)
        self.menubar.add_cascade(label="Help", menu=self.helpmenu)

        master.config(menu=self.menubar)

        # set control, graph and caption frames

        #self.framepreviewlabels = Frame(master)
        self.framepreviewcontrol = Frame(master)
        self.framepreviewcmd = Frame(master)
        self.frameplotcontrols = Frame(master)
        self.framecaption = Frame(master)
        self.framegraph = Frame(master)


        # choose variant

        v2=Frame(self.framepreviewcontrol, borderwidth=0)
        Label(v2, width=10, text = "RT parm" ).pack()
        self.selectvariant = Spinner.Spinner(v2,
                                        values=sorted(self.geninput.rtRange.keys(), reverse=True),
                                        width=10, init='TCLOUD', command=self.SetupPreview)
        self.selectvariant.pack()
        v2.pack(side='left')


        # choose the variant

        v3=Frame(self.framepreviewcontrol, borderwidth=0)
        Label(v3, text = "Description",  width=30).pack()
        self.description = Entry(v3, width=30)
        self.description.config(justify='center')
        self.description.pack()
        v3.pack(side='left')

        # set the range boxes

        v4=Frame(self.framepreviewcontrol, borderwidth=0)
        Label(v4, text = "Start", width=8).pack()
        self.rangestart = self.MakeEntryWidget(v4, width=8, Return=self.ShowPreview, MouseWheel=self.FixRange)
        self.rangestart.pack()
        v4.pack(side='left')

        v5=Frame(self.framepreviewcontrol, borderwidth=0)
        Label(v5, text="End", width=8).pack()
        self.rangefinis = self.MakeEntryWidget(v5, width=8, Return=self.ShowPreview, MouseWheel=self.FixRange)
        self.rangefinis.pack()
        v5.pack(side='left')

        self.SetRangeBoxes()

        # set number of samples

        v6=Frame(self.framepreviewcontrol, borderwidth=0)
        Label(v6, text="Num", width=4).pack()
        self.rangesamples = self.MakeEntryWidget(v6, init=11, width=4, Return=self.ShowPreview, MouseWheel=self.FixRange)
        self.rangesamples.pack()
        v6.pack(side='left')
        #print "rangesamples = ", self.rangesamples.get()

        # set the skewness

        v7=Frame(self.framepreviewcontrol, borderwidth=0)
        Label(v7, text="Skew", width=4).pack()
        self.rangeskew = self.MakeEntryWidget(v7, init=1, width=4, Return=self.ShowPreview, MouseWheel=self.FixRange)
        self.rangeskew.pack()
        v7.pack(side='left')
        # add the rt variant to the command file

        v9=Frame(self.framepreviewcontrol, width=15).pack(side='left')

        v8=Frame(self.framepreviewcontrol, borderwidth=3)
        self.runbtn = Button(v8, text='RUN', width=10,command=self.RunSBDART, borderwidth=2)
        self.runbtn.configure(background='pink')
        self.runbtn.pack(pady=5)
        self.addvariant = Button(v8, width=10, text='Add', command=self.AddRtParm)
        self.addvariant.configure(background='pink')
        self.addvariant.bind(self.rightmouse, self.NewRtParm)
        self.addvariant.pack(pady=5)
        v8.pack(side='left')

        self.abortbtn = Button(v8, width=10, text='Abort', command=self.AbortRun)

        # create command preview

        self.previewbox = Entry(self.framepreviewcmd)
        self.ShowPreview()
        self.previewbox.pack(fill='x', side='bottom', expand=1)

        # show legend checkbox

        self.showlegendVar = BooleanVar(value=True)
        self.showlegend = Checkbutton(self.frameplotcontrols, text="Show Legend", variable=self.showlegendVar,command=self.Plotit)
        self.showlegend.pack(side='left')

        # y-axis autoscale checkbox

        self.yautoscale = BooleanVar()
        self.yaxisautoscale = Checkbutton(self.frameplotcontrols, text="Y-Autoscale", variable=self.yautoscale,command=self.Plotit)
        self.yaxisautoscale.pack(side='left')

        # set up canvas

        self.fig = Figure(figsize=(9, 6), dpi=100) # 9,7
        self.fig.patch.set_facecolor('white')
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.framegraph)
        #self.canvas.mpl_connect('key_press_event', self.KeyPressHandler)
        #self.canvas.mpl_connect('button_press_event', self.ShowPlotData)
        self.canvas.mpl_connect('pick_event', self.SelectLine)
        self.canvas.mpl_connect('scroll_event', self.GraphMouseWheel)
        self.canvas.mpl_connect('motion_notify_event', self.HilightLine)

        self.caption = Text(self.framecaption, height=15, highlightbackground='black',borderwidth=1)
        self.caption.bind("<Double-Button-1>", self.SelectCaptionParm)
        self.caption.pack(fill='both', expand=1)

        self.framepreviewcontrol.pack(ipady=5)
        self.framepreviewcmd.pack(fill='x')
        self.frameplotcontrols.pack()
        self.framegraph.pack(expand=1)
        self.framecaption.pack(fill='both',expand=1)
        self.canvas._tkcanvas.pack(expand=1)

        # selectvariant here to ensure existence of rangeskew

        self.selectvariant.delete(0,"end")
        self.selectvariant.insert(0,"TCLOUD")
        self.SetupPreview()

        # load last file as a default

        try:
            with open("RecentFile", 'r') as myfile:
                filename = myfile.read()
            self.LoadFile(filename=filename)
            if not os.path.isfile(sbdartexe):
                self.Popup(self.framegraph, "Path to SBDART executable is not correct\nReplace line 25 of RunRT with correct path", wait=True)
        except:
            pass

    def SelectLine(self, event):
        '''show label of picked line in preview textbox'''
        lineinfo=""
        if isinstance(event.artist, Line2D):
            text = event.artist.get_label()
            if event.mouseevent.button == 1:
                txt = text.split(None, 1)[1]
                self.diffbase = txt
                self.PreviewLine("Baseline set to "+txt)
                self.Plotit()
            elif event.mouseevent.button == 2:
                for curve in self.ax.get_lines():
                    if curve.get_label() == text:
                        x=np.array(curve.get_xdata())
                        y=np.array(curve.get_ydata())
                        if self.parser.IOUT == 11:
                            xmin = x.min()
                            xmax = x.max()
                            quad = np.trapz(x, -y)
                            self.PreviewLine('{}:  Xmin={:.5g}  Xmax={:.5g}   Integral={:.5g}'.format(text,xmin,xmax,quad))
                        else:
                            ymin = y.min()
                            ymax = y.max()
                            quad = np.trapz(y, x)
                            self.PreviewLine('{}:  Ymin={:.5g}  Ymax={:.5g}   Integral={:.5g}'.format(text,ymin,ymax,quad))

                        break

    def GetPlotData(self, choice):
        ylbls0 = []
        ylbls1 = []
        ydata = []
        xdata = None
        title = self.ax.get_title()
        for curve in self.ax.get_lines():
            ylbls0.append(curve.get_label().split()[0])
            ylbls1.append(curve.get_label().split(' ',1)[1])
            xdata = curve.get_xdata()
            ydata.append(curve.get_ydata())
        if choice == 'View':
            w=11                # width of each column
            ws=str(w)
            nc = len(ylbls0)
            tpad = (w*(nc+1)-len(title))/2
            txt = ' '*tpad + title +'\n\n'
            fm1 = ' ' * w + ('\t{:' + ws + '}') * nc + '\n'
            fm2 = '\n{:'+ws+'.3e}' + ('\t{:'+ws+'.3e}') * nc
            txt += fm1.format(*ylbls0)
            txt += fm1.format(*ylbls1)
            for i, x in enumerate(xdata):
                y = []
                for j in range(0, len(ydata)):
                    y.append(float(ydata[j][i]))
                txt += fm2.format(x, *y)
            return txt
        elif choice == 'Copy':
            txt = 'x=[' + ','.join(['{:11.3e}'.format(float(x)) for x in xdata]) + ']\n'
            yparts = []
            for yy in ydata:
                t = '[' + ','.join(['{:11.3e}'.format(float(y)) for y in yy]) + ']'
                yparts.append(t)
            txt += 'y=[' + ',\n'.join(yparts) + ']'
            return txt

    def PlotData(self, choice):
        if choice == 'Save':
            for k in range(0, 100):
                dir = ""
                if os.path.isdir("RUNS"):
                    dir = "RUNS{}".format(os.path.sep)
                plotname = "{}{}_{}.png".format(dir, self.runname, str(k).zfill(2))
                if not os.path.isfile(plotname):
                    self.fig.savefig(plotname)
                    self.Popup(self.framegraph,"{} saved".format(plotname))
                    return
            self.Popup(self.framegraph, "Could not write " + plotname)

        else:
            txt = self.GetPlotData(choice)
            if choice == 'View':
                self.Popup(self.framegraph,txt)
            r = Tk()
            r.withdraw()
            r.clipboard_clear()
            r.clipboard_append(txt)
            r.destroy()

    def HilightLine(self, event):
        '''hilight line on mouse hover'''
        hilightwidth = 2
        normalwidth = 1
        hotone = ''
        if not self.ax.lines:
            return
        for curve in self.ax.get_lines():
            if curve.contains(event)[0] and not hotone:
                hotone = curve.get_label()
                curve.set_linewidth(hilightwidth)
            else:
                curve.set_linewidth(normalwidth)

        for legline in self.leg.get_lines():
            if legline.get_label() == hotone:
                legline.set_linewidth(hilightwidth)
            else:
                legline.set_linewidth(normalwidth)

        self.fig.canvas.draw()

    def SelectCaptionParm(self, event):
        '''
        select a parameter in the caption box with a double click
        show that parameter in the selectvariant entry box
        :param event:
        :return:
        '''
        pnt = self.caption.index("@%d,%d" % (event.x, event.y))
        lind = pnt.split('.')[0]
        pnt0="{}.0".format(lind)
        pnt1="{}.0-1c".format(int(lind)+1)
        line = self.caption.get(pnt0,pnt1)
        #print "line ",line
        if line.count('='):
            parm,vals=line.split('=')
            self.selectvariant.delete(0,"end")
            self.selectvariant.insert(0,parm.upper().strip())
            self.SetupPreview()

    def ZoomAxis(self):
        '''
        adjust x and y axis limits
        :return:
        '''
        xmn,xmx=self.plotwindow[0:3:2]
        xtic=self.NiceStep(xmx-xmn)
        xlim=[xtic*round(xmn/xtic), xtic*round(xmx/xtic)]
        self.ax.set_xlim(xlim)
        ymn,ymx=self.plotwindow[1:4:2]
        ytic=self.NiceStep(ymx-ymn)
        ylim=[ytic*round(ymn/ytic), ytic*round(ymx/ytic)]
        self.ax.set_ylim(ylim)

    def GraphMouseWheel(self, event):
        """
        interprets mouse roll events to change axis limits
        rolls below the x-axis controls the x-axis
        rolls left of the y-axis controls the y-axis
        """
        if self.plottype == 'xy':
            self.GraphMouseWheelXY(event)
        else:
            self.GraphMouseWheelPolar(event)

    def GraphMouseWheelXY(self, event):
        x=event.x
        y=event.y
        (xbmn,ybmn,xbmx,ybmx)=self.rectXY
        xsz,ysz=self.canvas.get_width_height()
        xpn=float(x)/xsz
        ypn=float(y)/ysz

        if (xpn-xbmn)*(ypn-ybmn) > 0:
            self.plotwindow=[]
            self.Plotit()
            return

        if self.plotwindow:
            xmn,ymn,xmx,ymx=self.plotwindow
        else:
            xmn,xmx = self.ax.get_xlim()
            ymn,ymx = self.ax.get_ylim()
        xtic=0.1*(xmx-xmn)
        ytic=0.1*(ymx-ymn)
        ds = 1 if event.step>0 else -1
        if ypn < ybmn:                      # adjust x-axis
            xnorm = (xpn-xbmn)/(xbmx-xbmn)
            if xnorm < 0.5:
                xmn += ds*xtic
            else:
                xmx += ds*xtic
        elif xpn < xbmn:                    # adjust y-axis
            ynorm = (ypn-ybmn)/(ybmx-ybmn)
            if ynorm < 0.5:
                ymn += ds*ytic
            else:
                ymx += ds*ytic
        self.plotwindow=[xmn,ymn,xmx,ymx]
        self.Plotit()

    def GraphMouseWheelPolar(self, event):
        x=event.x
        y=event.y
        (xbmn,ybmn,xbmx,ybmx)=(0.125,0.1,0.6,0.9)  # need a programatic way to get bounding box in normed coords
        xsz,ysz=self.canvas.get_width_height()
        xpn=float(x)/xsz
        ypn=float(y)/ysz

        if xpn < xbmx:
            self.colorbarzoom = []
            self.Plotit()
        else:
            vmn,vmx=(0,0)
            if self.colorbarzoom:
                vmn,vmx= self.colorbarzoom
                print 'vmn vmx {}  {}'.format(vmn,vmx)
                ztic = 0.1*(vmx-vmn)
                if ztic > 0:
                    ynorm = (float(y)/ysz-ybmn)/(ybmx-ybmn)
                    ds = 1.0 if event.step > 0 else -1.0
                    if ynorm < 0.5:
                        vmn += ds * ztic
                    else:
                        vmx += ds * ztic

            self.colorbarzoom = [vmn,vmx]
            self.Plotit()

    def MakeEntryWidget(self, frame, **kwargs):
        """
        create an entry widget with standard callbacks references
        :param self:
        :param frame:  parent frame
        :param kwargs:
            init        -- initial value
            width       -- width of entry widget
            Return      -- call back function bound to return key
            MouseWheel  -- call back function bound to mousewheel up or down
            side        -- side to pack
        :return:
        """

        if kwargs.has_key('width'):
            entry=Entry(frame, width=kwargs['width'])
        else:
            entry=Entry(frame)

        for k in kwargs.keys():
            #print "MakeEntry {} : {}".format(k, kwargs[k])
            if k == 'init':
                entry.insert(0,kwargs[k])
            if k == 'Return':
                entry.bind('<Return>', kwargs['Return'])
            if k == 'MouseWheel':
                entry.bind('<MouseWheel>', kwargs['MouseWheel'])

        return entry

    def HelpGUI(self):
        if not self.rtdocwindow:
            fh = open("runrtdoc.txt", "r")
            lines = fh.readlines()
            msg = "".join(lines)
            fh.close()
            self.Popup(self.framegraph, msg)

    def ViewOutput(self):
        msg=self.sbdartoutput
        if msg:
            self.Popup(self.framegraph, msg)

    def HelpTopic(self):
        key = self.selectvariant.get()
        pattern = " " + key + ":"
        try:
            txt = self.rtdoctext
            ind = txt.search(pattern, "1.0", stopindex="end")
            txt.mark_set("insert", ind)
            txt.see("insert")
        except:
            self.rtdocwindow = None

    def HelpRtDoc(self):
        if not self.rtdocwindow:
            fh = open("rtdoc.txt", "r")
            lines = fh.readlines()
            msg = "".join(lines)
            fh.close()
            self.rtdoctext, self.rtdocwindow = self.Popup(self.framegraph, msg)

    def ShowFilePopup(self, file):
        if not self.rtdocwindow:
            fh = open(file, "r")
            lines = fh.readlines()
            msg = "".join(lines)
            fh.close()
            self.rtdoctext, self.rtdocwindow = self.Popup(self.framegraph, msg)

    def ShowChoice(self, event):
        var = event.widget.get().split('=')[0]
        cycleseq = self.geninput.CycleSequence()
        self.variant_to_plot = cycleseq.index(var)
        self.SetGroupMenu()
        self.PreviewLine("")
        try:
            self.Plotit()
        except:
            print "Error in ShowChoice"
            print "var:             ", var
            print "variant_to_plot: ", self.variant_to_plot
            print "cycleseq:        ", cycleseq

    def SetRtMenu(self):
        """
        replace rtmenu with parms appropriate to sbdart run
        """
        menu=self.fmenu
        parser=self.parser
        menu.delete(0, 'end')
        seperator = True
        for k in sorted(parser.menucheck.keys()):
            if parser.rtdefault == k:
                parser.menucheck[k].set(True)
            if seperator and not k.startswith(' '):
                menu.add_separator()
                seperator = False
            menu.add_checkbutton(label=k, variable=parser.menucheck[k], onvalue=True, offvalue=False, command=lambda key=k: self.PlotGroup(key))

    def SetGroupMenu(self):
        """
        Create Plot Group menu items
        suppose cmd sequence is
                         A=0;1
                         B=2;3
                         C=4;5
                         D=6;7

        iout=20,21:     all variations show up in groupmenuSpinners
        iout=10:        a bunch of chkboxes are created for all C values
                        two spinboxes are created to select the B and D values
        iout=1,2,11:    a bunch of chkboxes are created for all B values
                        two spinboxes are created to select the "C" and "D" values
        :return:
        """
        # self.plot_sequence = tags
        # consider tags as:
        # tags=[['A=0','A=1'],['B=2','B=3], ['C=4','C=5'], ['D=6', 'D=7']
        ichk = self.variant_to_plot
        tags = self.geninput.GetParmMenuItems()
        if self.parser.IOUT == 10:
            tags=tags[1:]                # tags=[['B=2','B=3], ['C=4','C=5'], ['D=6', 'D=7']
        elif self.parser.IOUT in (20,21):
            ichk=-1                      # all variations show up in groupmenuSpinners
        self.groupmenu.delete(0, 'end')
        # print "iout:     ",self.parser.IOUT
        # print "ichk:     ",ichk
        # print "tags:     ",tags

        self.groupmenuChkBox = OrderedDict()
        if len(tags) > 0 and not ichk == -1:
            for choice in tags[ichk]:  #
                self.groupmenuChkBox[choice] = BooleanVar(value=True)
                cleaned_choice = self.GetLinelabel(choice)
                self.groupmenu.add_checkbutton(label=cleaned_choice, variable=self.groupmenuChkBox[choice], onvalue=True,
                                               offvalue=False, command=self.Plotit)
        for spinbox in self.groupmenuSpinners:
            spinbox.destroy()

        self.groupmenuSpinners = []
        jchk = -1
        rbutn = "<Button-2>" if platform.system() == "Darwin" else "<Button-3>"  # OSX version reverses Button-2 and Button-3

        for choices in tags:  # choices = ['B=0.0','B=10', 'B=20']
            jchk += 1
            if jchk == ichk: continue
            spinbox = Spinbox(self.frameplotcontrols, width=10, values=choices, command=self.Plotit, repeatdelay=500)
            spinbox.pack(side='left')
            spinbox.bind(rbutn, self.ShowChoice)
            spinbox.bind('<MouseWheel>', self.SpinBoxIncrement)
            self.groupmenuSpinners.append(spinbox)

    def SetupNominalPlots(self):
        for spinbox in self.optionSpinners:
            spinbox.destroy()
        self.optionSpinners = []
        self.Plotit()

    def ValidForEphemeris(self, constantParms, parms):
        '''
        check if run parameters are consistent with ephemeris plots
        :param constantParms:
        :param parms:
        :return:
        '''
        iout = 10
        for c in constantParms:
            c = c.split('#')[0]
            if c.startswith('IOUT'):
                iout = int(c.split('=')[1])
        msg = ''
        if iout != 10:
            msg = "* This plot option requires output generated with IOUT=10\n\n"
        if not (parms.startswith('SZA') or parms.startswith('CSZA')):
            msg += "* Diurnal variation plots require that SZA is the first variable parameter\n\n"

        if msg:
            self.Popup(self.framegraph, msg)
            self.optionDiurnalPlot.set(0)
            return False
        else:
            return True

    def SetupHourly(self):
        '''
        generate spinbox input widgets for day, lat used with solar ephemeris plots
        :return:
        '''

        variableParms, constantParms, parms = self.GetParmsInCmd(clean=True)

        valid = self.ValidForEphemeris(constantParms, parms)

        if valid:

            for spinbox in self.optionSpinners:
                spinbox.destroy()
            self.optionSpinners=[]

            spinbox = self.SpinboxEphem('day', range(1,366,5), 171)
            self.optionSpinners.append(spinbox)
            spinbox = self.SpinboxEphem('lat', range(-90,91,5), 35)
            self.optionSpinners.append(spinbox)

            self.Plotit()

    def SpinboxEphem(self, parm, span, default, **kwargs):
        '''

        :param parm:     variable name string (e.g., 'day')
        :param span:     range of values      (e.g., range(1,366,5)
        :param default:  starting value       (e.g., 171
        :param kwargs:
        :return:        spinbox object
        '''
        kwargs['values'] = ['{}={}'.format(parm, v) for v in span]
        kwargs['repeatdelay']=500
        kwargs['wrap']=True
        kwargs['command']=self.Plotit
        wid = 0
        for v in kwargs['values']:
            w = len(v)
            if w > wid:
                wid = w
        kwargs['width'] = wid
        spinbox = Spinbox(self.frameplotcontrols, **kwargs)
        if default:
            spinbox.delete(0,'end')
            spinbox.insert(0, '{}={}'.format(parm, default))
        spinbox.pack(side='left')
        spinbox.bind('<MouseWheel>', self.SpinBoxIncrement)
        return spinbox

    def SetupHourlyPlace(self):
        '''
        generate spinbox input widgets for day, lat and lon used with solar ephemeris plots
        :return:
        '''

        variableParms, constantParms, parms = self.GetParmsInCmd(clean=True)

        valid = self.ValidForEphemeris(constantParms, parms)

        if valid:

            for spinbox in self.optionSpinners:
                spinbox.destroy()
            self.optionSpinners=[]

            spinbox = self.SpinboxEphem('day', range(1,366), 171)
            self.optionSpinners.append(spinbox)
            spinbox = self.SpinboxEphem('lat', range(-90,91), 35)
            self.optionSpinners.append(spinbox)
            spinbox = self.SpinboxEphem('lon', range(-180,181), -120)
            self.optionSpinners.append(spinbox)

            self.Plotit()

    def SetupDailySpinner(self, parm):
        variableParms, constantParms, parms = self.GetParmsInCmd(clean=True)
        valid = self.ValidForEphemeris(constantParms, parms)

        if valid:

            for spinbox in self.optionSpinners:
                spinbox.destroy()
            self.optionSpinners=[]

            if parm == 'day':
                spinbox = self.SpinboxEphem('day', range(1,366,5), 177)
            else:
                spinbox = self.SpinboxEphem('lat', range(-90,91,5), 35)
            self.optionSpinners.append(spinbox)
            self.Plotit()


    def GetOptionSpinnerValues(self):
        values = []
        for spinbox in self.optionSpinners:
            values.append(spinbox.get())
        return values

    def SpinBoxIncrement(self, event):
        """
        increment Spinbox selection using mousewheel motion
        :param event:
        :return:
        """
        if event.delta > 0:
            event.widget.invoke('buttonup')
        elif event.delta < 0:
            event.widget.invoke('buttondown')

    def GetPlotSequence(self):
        """
        get a sequence of plot keys based on groupmenu and fmenu selections
        the plot key is the dictionary key into yvariable
        cklist is the list of second order variants
        pklist is the list of third and higher order variants
        The method used to construct the yvariable key must be consistent with RtReader.keymaker
        :return: dictionary of activated plot quantities. Value entries are [marker, color, linelabel, plottitle]
        """
        # construct keys from checked elements of groupmenuChkBox
        markers = ',.ov^<>Dd1234sp*hH+x|_'
        colors = 'bgrcmyk'
        constr = ""
        cklist = []
        try:
            for key in self.groupmenuChkBox.keys():
                if self.groupmenuChkBox[key].get() and not key.endswith('#'):
                    cklist.append(key)
        except:
            pass

        pklist = []
        try:
            for sb in self.groupmenuSpinners:
                seq = sb.get()
                pklist.append(seq)
            constr = " ".join(sorted(pklist))
        except:
            constr = ""

        # print "GetPlotSequence constr=", constr
        #
        # print "GetPlotSequence cklist=", cklist

        fkeys = self.parser.Plotables()

        pdict = OrderedDict()
        j = -1
        for fkey in fkeys:
            j = (j + 1) % len(markers)  # markers used for flux types
            if self.parser.IOUT in (20,21):
                key = "{} {}".format(fkey, constr) if constr else fkey
                linelabel = key
                plottitle = constr
                pdict[key] = [colors[0], markers[0], linelabel, plottitle]
            else:
                if len(cklist) == 0:
                    pdict[fkey] = [colors[0], markers[j], fkey, ""]
                else:
                    p = -1
                    for pkey in cklist:
                        if not self.groupmenuChkBox[pkey].get(): continue
                        p = (p + 1) % len(colors)  # colors for different parm selections
                        rkey = "{} {}".format(pkey, constr) if constr else "{}".format(pkey, fkey)
                        rkey = " ".join(sorted(rkey.split()))  # sorted => keys independent of nesting order
                        key = "{} {}".format(fkey, rkey)
                        linelabel = "{} {}".format(fkey, pkey)
                        plottitle = constr
                        pdict[key] = [colors[p], markers[j], linelabel, plottitle]

        return pdict

    def SetRange(self, vec, oldrange):
        minv = min(vec + oldrange[0:1])
        maxv = max(vec)
        if maxv > 0:
            maxv *= 1.1
        maxv = max([maxv, oldrange[1]])
        return [minv, maxv]

    def MakeAx(self,ptype):
        """
        change to given plot type, either xy or polar, if necessary
        :param ptype:
        :return:
        """
        if ptype == self.plottype and ptype == 'xy':
            return
        self.plottype = ptype
        self.fig.clear()
        if ptype == "polar":
            self.ax = self.fig.add_subplot(111, projection='polar')
            self.fig.tight_layout(rect=self.rectPolar)

        else:
            self.ax = self.fig.add_subplot(111)
            self.fig.tight_layout(rect=self.rectXY)

    def PrintArray(self,array):
        nx,ny=array.shape
        for iy in range(0,ny):
            line="".join(["{:10.3f}".format(x) for x in array[:,iy]])
            print line

    def PlotGroup(self, key):
        """
        performs RT menu book-keeping.  Switches off checkbox entries with conflicting units
        call Plotit after book-keeping
        :param key:
        :return:
        """
        if not key.startswith(' '):
            if self.parser.IOUT in (20,21):
                for k in self.parser.menucheck.keys():
                    if k.startswith(' '):
                        continue
                    if k != key:
                        self.parser.menucheck[k].set(False)
            else:
                units = self.parser.rtunits[key]
                for k in self.parser.menucheck.keys():
                    if k.startswith(' '):
                        continue
                    if units != self.parser.rtunits[k]:
                        self.parser.menucheck[k].set(False)
        self.Plotit()

    def GetLinelabel(self, label):
        '''
        remove psuedo variants from line label, replace underscores with white space
        :param label:
        :return: filtered line label
        '''
        s = []
        for lbl in label.split():
            if not lbl[0].isalpha():
                lbl = lbl.split('=')[1].replace('_',' ')
            s.append(lbl)
        return " ".join(s)

    def goodkey(self, key, dic=None):
        diction = dic if dic else self.yvariable
        if key in diction.keys():
            return True
        else:
            frame, filename,line_number, function_name, lines, index = inspect.getouterframes(inspect.currentframe())[1]
            msg = "A bad dictionary key was found at line number "
            msg += "{} of {} in {}\n\n".format(line_number, os.path.basename(filename), function_name)
            msg += "\n".join(lines)
            msg += "\nkey='{}' is not a member of \n\n".format(key)
            for ky in diction.keys():
                msg+="{}\n".format(ky)
            self.Popup(self.framegraph, msg)
            return False

    def Plotit(self):

        if not self.diffbase and self.optionComparisonPlot.get() > 0:
            self.optionComparisonPlot.set(0)
            self.Popup(self.framegraph,"No comparison base set")
            return

        if self.parser.IOUT in (20,21):
            self.Plotit20()
        elif self.parser.IOUT == 11:
            self.Plotit11()
        elif self.parser.IOUT in (1,2):
            self.Plotit01()
        else:
            self.Plotit10()

    def Plotit01(self):
        self.ClearPlot()
        self.MakeAx('xy')
        pdict = self.GetPlotSequence()

        # for ykey in self.yvariable.keys():
        #     print "ykey: ",ykey
        #     print self.yvariable[ykey]

        color =''
        plottitle=''
        basequant = self.diffbase if self.optionComparisonPlot.get() else ''
        for pkey in pdict.keys():
            color = pdict[pkey][0]
            marker = pdict[pkey][1]
            linelabel = self.GetLinelabel(pdict[pkey][2])
            plottitle = pdict[pkey][3]
            rtkey=pkey.split()[0]
            wlkey = pkey.replace(rtkey, "WL", 1)
            if not wlkey in self.yvariable.keys():
                wlkey = "WL"
            ylabel = self.parser.rtunits[rtkey]
            xlabel = "$Wavelength (\mu m)$"
            iswavenumber = self.parser.menucheck.has_key(' Wavenumber') and self.parser.menucheck[' Wavenumber'].get()
            isefftemp = self.parser.menucheck.has_key(' EffectiveTemp') and self.parser.menucheck[' EffectiveTemp'].get()
            istransmit = self.parser.menucheck.has_key(' Transmitted') and self.parser.menucheck[' Transmitted'].get()
            x = self.yvariable[wlkey][:]
            self.xvariable = x
            if self.goodkey(pkey):
                y = self.yvariable[pkey][:]

                if basequant:
                    diffquant = linelabel[linelabel.find(' ')+1:]
                    if basequant == diffquant: continue
                    diffquant=diffquant.replace(' ','_')
                    basequant=basequant.replace(' ','_')
                    bkey = pkey.replace(diffquant, basequant)
                    plottitle = "Difference from {}  {}".format(basequant, plottitle)
                    if self.goodkey(bkey):
                        if self.optionComparisonPlot.get() == 1:
                            y = [a-b for a,b in zip(y,self.yvariable[bkey][:])]
                        elif self.optionComparisonPlot.get() == 2:
                            y = [a/b if b>0 else 1e-6 for a,b in zip(y,self.yvariable[bkey][:])]
                            ylabel = 'Ratio'

                if istransmit:
                    tpkey = pkey.replace(rtkey, "TOPDN")
                    ylabel=''
                    linelabel=linelabel.replace(rtkey, '{}/{}'.format(rtkey, 'TOPDN'))
                    y = [yy / yyt if yyt > 0 else 0 for yy,yyt in zip(y, self.yvariable[tpkey][:])]

            else:
                return

            if isefftemp:
                ylabel = "Effective Temperature (K)"
                for i in range(0, len(y)):
                    y[i] = self.PlanckTemp(x[i],y[i])

            if iswavenumber:
                xlabel="$Wavenumber (cm^{-1})$"
                if ylabel:
                    ylabel=ylabel.replace("\mu m","cm^{-1}") if not isefftemp else "Effective Temperature (K)"
                x = [1e4/xx for xx in x]
                if not isefftemp and ylabel:
                    for i in range(0, len(y)):
                        y[i] = 1e4*y[i]/x[i]**2
            self.ax.plot(x, y, color=color, marker=marker, label=linelabel, markersize=3, picker=5)
            self.yrange = self.SetRange(y, self.yrange)

        if not color:
            return
        if self.yautoscale.get():
            self.plotwindow=[]
            self.yrange = [float('inf'), float('-inf')]
        else:
            if self.plotwindow:
                self.ZoomAxis()
            else:
                self.ax.set_ylim(self.yrange)

        self.ShowLegend(len(pdict.keys()))
        self.ax.set_title(plottitle)
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        self.canvas.show()

    def PlanckTemp(self, w, r):
        '''
        compute the effective temperature given wavelength and irradiance
        B = 2hc^2/w^5  (1/e^(hc/wkt)-1)  (w/sr/m^3)
        the irradiance Bi=pi * B
        (2*pi*h*c^2)/(Bi*w^5)+1 = exp(hc /w k t)
        t=hc/(wk*log(lhs))
        :param w:  wavelength (um)
        :param r:  irradiance (w/m2/um)
        :return: effective temperature (K)
        '''
        wmks = w*1e-6
        Bi = r*1e6
        h=6.626e-34
        c=2.9979e8
        k=1.38065e-23
        lhs = 1+2*np.pi*h*c**2/(Bi*wmks**5)
        efftemp = h*c/(wmks*k*np.log(lhs))
        return efftemp


    def Plotit10(self):
        self.ClearPlot()
        self.MakeAx('xy')
        pdict = self.GetPlotSequence()
        color = ""
        plottitle = ""
        xlabel = "$"+self.xlabel+"$"
        x = np.array([float(xx) for xx in self.xvariable])
        diurnalop = self.optionDiurnalPlot.get()

        if diurnalop in [1, 2]:
            day, lat, lon = self.GetEphemInputs()
            date, noon, daylight, solfac, wsza, x = self.SetupEphemHours(x, day, lat, lon)
            wfloor = np.floor(wsza)
            wt = wsza - wfloor
            isza = wfloor.astype(int)
            iszap = np.ceil(wsza).astype(int)
            xlabel = 'UTC TIME (hours)'
            noontime = "{}:{:02}:{:02}".format(int(noon), int(noon*60) % 60, int(noon*3600) % 60)
            plottitle = 'Date {}   noon={}  daylight={:.2f}'.format(date, noontime, daylight)

        elif diurnalop == 3:
            xlabel = 'Latitude (deg)'
            day, junk, junk = self.GetEphemInputs()
            lats = np.arange(-90,91,5)
            date, solfacarr, timearr, wszaarr = self.DiurnalAverageSetup(day, lats, x)
            plottitle = 'Date {}'.format(date)
        elif diurnalop == 4:
            xlabel = 'Day Number'
            junk, lat, junk = self.GetEphemInputs()
            days = np.arange(1,366,5)
            date, solfacarr, timearr, wszaarr = self.DiurnalAverageSetup(days, lat, x)

        basequant = self.diffbase if self.optionComparisonPlot.get() else ''
        for pkey in pdict.keys():
            color = pdict[pkey][0]
            marker = pdict[pkey][1]
            linelabel = self.GetLinelabel(pdict[pkey][2])
            if diurnalop == 0:
                plottitle = pdict[pkey][3]
            rtkey=pkey.split()[0]
            ylabel = self.parser.rtunits[rtkey]
            intensity = self.parser.menucheck[" Intensity"].get()
            transmitted = self.parser.menucheck[" Transmitted"].get()

            if self.goodkey(pkey):
                y = np.array(self.yvariable[pkey][:])
                if transmitted:
                    intensity=False
                    ylabel=''
                    topkey = pkey.replace(rtkey,'TOPDN')
                    linelabel=linelabel.replace(rtkey, '{}/{}'.format(rtkey, 'TOPDN'))
                    yt = np.array(self.yvariable[topkey][:])
                    y = np.where(yt>0, y/yt, np.zeros(len(y)))

                if diurnalop in [1, 2]:
                    y = solfac*(y[isza]*(1-wt)+y[iszap]*wt)
                elif diurnalop == 3:
                    x, y = self.DiurnalAverage(lats, solfacarr, timearr, wszaarr, y)
                elif diurnalop == 4:
                    x, y = self.DiurnalAverage(days, solfacarr, timearr, wszaarr, y)

                if basequant:
                    diffquant = linelabel[linelabel.find(' ')+1:]
                    if basequant == diffquant: continue
                    if diurnalop == 0:
                        plottitle = "Difference from {}    {}".format(basequant, plottitle)
                    diffquant=diffquant.replace(' ','_')
                    basequant=basequant.replace(' ','_')
                    bkey = pkey.replace(diffquant, basequant)
                    if self.goodkey(bkey):
                        yb = np.array(self.yvariable[bkey])
                        if transmitted:
                            topkey = bkey.replace(rtkey,'TOPDN')
                            #print 'bkey={}  topkey={}'.format(bkey, topkey)
                            ytb = np.array(self.yvariable[topkey][:])
                            yb = np.where(ytb>0, yb/ytb, np.zeros(len(yb)))

                        if diurnalop in [1, 2]:
                            yb = solfac*(yb[isza]*(1-wt)+yb[iszap]*wt)
                        elif diurnalop == 3:
                            junk, yb = self.DiurnalAverage(lats, solfacarr, timearr, wszaarr, yb)
                        elif diurnalop == 4:
                            junk, yb = self.DiurnalAverage(days, solfacarr, timearr, wszaarr, yb)

                        if self.optionComparisonPlot.get() == 1:
                            y -= yb
                        elif self.optionComparisonPlot.get() == 2:
                            y /= np.where(yb>0, yb, 1e-6*y)
                            ylabel = 'Ratio'

                if intensity:
                    ylabel=self.parser.IntensityLabel(ylabel)
                    ewkey = pkey.replace(rtkey, "FFEW")
                    ffew = self.yvariable[ewkey][0]
                    y /= ffew

                self.ax.plot(x, y, color=color, marker=marker, label=linelabel, markersize=3, picker=5)

                self.yrange = self.SetRange([y.min(),y.max()], self.yrange)
            else:
                return

        if not color:
            return
        if self.yautoscale.get():
            self.plotwindow=[]
            self.yrange = [float('inf'), float('-inf')]
        else:
            if self.plotwindow:
                self.ZoomAxis()
            else:
                self.ax.set_ylim(self.yrange)

        self.ShowLegend(len(pdict.keys()))
        self.ax.set_title(plottitle)
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        self.canvas.show()

    def DiurnalAverageSetup(self, days, lats, x):
        wszaarr = []
        solfacarr = []
        timearr = []

        if hasattr(lats,"__len__"):
            days = [days]*len(lats)
        else:
            lats = [lats]*len(days)
        for k in range(len(days)):
            day = days[k]
            lat = lats[k]
            date, noon, daylight, solfac, wsza, time = self.SetupEphemHours(x, day, lat, 0)
            wszaarr.append(wsza)
            solfacarr.append(solfac)
            timearr.append(time)
        return date, solfacarr, timearr, wszaarr

    def DiurnalAverage(self, lats, solfacarr, timearr, wszaarr, y):
        '''
        compute diurnal average of irradiance quantity y
        :param lats:  independent variable, either lats or days
        :param solfac: solar orbital factor
        :param solfacarr:  array of solar orbital factors
        :param timearr:    array of times
        :param wszaarr:    array of floating point indicies
        :param x:
        :param y:
        :return:
        '''
        farray = []
        for k in range(len(lats)):
            wsza = wszaarr[k]
            wfloor = np.floor(wsza)
            wt = wsza - wfloor
            isza = wfloor.astype(int)
            iszap = np.ceil(wsza).astype(int)
            solfac = solfacarr[k]
            t = timearr[k]
            f = solfac * (y[isza] * (1 - wt) + y[iszap] * wt)
            farray.append(np.trapz(f, x=t) / 24)
        y = np.array(farray)
        x = np.array(lats)
        return x, y

    def Plotit11(self):
        self.ClearPlot()
        self.MakeAx('xy')
        pdict = self.GetPlotSequence()

        plottitle = ""

        if len(pdict.keys()) == 0:
            return

        basequant = self.diffbase if self.optionComparisonPlot.get() else ''
        for pkey in pdict.keys():
            color = pdict[pkey][0]
            marker = pdict[pkey][1]
            linelabel = self.GetLinelabel(pdict[pkey][2])
            plottitle = pdict[pkey][3]
            rtkey=pkey.split()[0]
            phikey = pkey.replace(rtkey, "FFEW")
            zkey=pkey.replace(rtkey, "ZZ")
            if not zkey in self.yvariable.keys():
                zkey = 'ZZ'
            y = self.yvariable[zkey]
            ylabel = "Altitude (km)"
            xlabel = self.parser.rtunits[rtkey]
            intensity = self.parser.menucheck[" Intensity"].get()
            if self.goodkey(pkey):
                x = self.yvariable[pkey][:]
                self.xvariable = x

                if basequant:
                    diffquant = linelabel[linelabel.find(' ')+1:]
                    if basequant == diffquant: continue
                    plottitle = "Difference from {}  {}".format(basequant, plottitle)
                    diffquant=diffquant.replace(' ','_')
                    basequant=basequant.replace(' ','_')
                    bkey = pkey.replace(diffquant,basequant)
                    if self.goodkey(bkey):
                        if self.optionComparisonPlot.get() == 1:
                            x = [a-b for a,b in zip(x,self.yvariable[bkey][:])]
                        elif self.optionComparisonPlot.get() == 2:
                            x = [a/b if b>0 else 1e-6 for a,b in zip(x,self.yvariable[bkey][:])]
                            xlabel = 'Ratio'

                    else:
                        return

                if intensity:
                    xlabel=self.parser.IntensityLabel(xlabel)
                    for i in range(0, len(y)):
                        ffew = self.yvariable[phikey][i]
                        x[i] /= ffew
                self.ax.plot(x, y, color=color, marker=marker, label=linelabel, markersize=3, picker=5)
                self.yrange = self.SetRange(x, self.yrange)
            else:
                return

        if self.yautoscale.get():
            self.plotwindow=[]
            self.yrange = [float('inf'), float('-inf')]
        else:
            if self.plotwindow:
                self.ZoomAxis()
            else:
                self.ax.set_xlim(self.yrange)

        self.ShowLegend(len(pdict.keys()))
        self.ax.set_title(plottitle)
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        self.canvas.show()

    def Plotit20(self):
        azimuths=np.radians(self.parser.phi)
        zeniths = np.array(self.parser.zen)
        theta,phi=np.meshgrid(zeniths, azimuths)
        pdict = self.GetPlotSequence()
        values=[]
        self.ClearPlot()
        intensity = self.parser.menucheck[" Intensity"].get()
        if self.yautoscale.get():
            self.colorbarzoom=[]

        for pkey in pdict.keys():
            klist = pkey.split()
            rtparms = " ".join(klist[1:])
            rtkey = klist[0]
            if self.goodkey(rtkey, self.parser.rtunits):
                radkey = "RADIANCE {}".format(rtparms)
                units=self.parser.rtunits[rtkey]
                if self.goodkey(radkey) and self.goodkey(pkey):
                    values=np.array(self.yvariable[radkey])
                    label="{}  ({})".format(rtkey.capitalize(),self.parser.rtunits[rtkey])
                    if rtkey != "RADIANCE":
                        values /= self.yvariable[pkey]
                        label="Radiance/{}  ({})".format(rtkey.capitalize(),self.parser.rtunits[rtkey])
                else:
                    return
                if intensity and rtkey == "RADIANCE":
                    label = self.parser.IntensityLabel(label)
                    ewkey = pkey.replace(rtkey, "FFEW")
                    ew=self.yvariable[ewkey]
                    values = values/ew
                self.MakeAx("polar")
                self.ax.set_theta_offset(np.pi/2)
                self.ax.set_theta_direction(-1)
                if self.colorbarzoom:
                    vmin, vmax = self.colorbarzoom
                    if vmin == 0 and vmax == 0:
                        self.colorbarzoom = [values.min(), values.max()]
                    else:
                        values = values.clip(min=vmin,max=vmax)
                cp=self.ax.contourf(phi, theta, values,20)
                self.fig.colorbar(cp, ax=self.ax, orientation='vertical',pad=0.1,label=label)
                #self.fig.tight_layout(rect=[0.01,0.01,0.9,0.99])
                self.canvas.show()
            else:
                return



    def SetupEphemHours(self, x, day, lat, lon):
        '''
        set up parameters for hourly ephemeris computation
        :param x:
        :return:
        date        date string
        noon        time at which sun reaches highest elevation at this day and location
        daylight    hours of daylight
        solfac      irradiance factor cos(sunzen)/rsun**2
        wsza        floating point index for yvariable interpolation
        hours       time array
        '''
        ephm = Ephemeris.SolarEphemeris(lat, lon)
        noon, daylight = ephm.suntimes(day)
        dayhours = max([daylight, 8])
        hours = np.linspace(noon - 0.5 * dayhours, noon + 0.5 * dayhours, 40)
        zen, phi, solfac = ephm.sunpos(day, hours)
        date = str(datetime.date.fromordinal(day))[5:]
        sza = x * np.pi / 180
        wsza = np.interp(zen, sza, np.arange(0, len(sza)), left=0, right=len(x) - 1)
        return date, noon, daylight, solfac, wsza, hours

    def GetEphemInputs(self):
        '''
        get the parameter values associated with solar ephemeris options
        :return:
        '''
        optionvalues = self.GetOptionSpinnerValues()
        lat = lon = day = 0
        for op in optionvalues:
            if op.startswith('lat'):
                lat = int(op.split('=')[1])
            if op.startswith('lon'):
                lon = int(op.split('=')[1])
            if op.startswith('day'):
                day = int(op.split('=')[1])
        return day, lat, lon

    def ShowLegend(self, numleg):
        if not self.showlegendVar.get():
            return
        if numleg < 5:
            self.leg = self.ax.legend(loc='best', fontsize='medium')
        elif numleg < 10:
            self.leg = self.ax.legend(loc='best', fontsize='small')
        elif numleg < 20:
            self.leg = self.ax.legend(loc='best', fontsize='x-small')
        else:
            self.leg = self.ax.legend(loc='best', fontsize='xx-small')

    def Popup(self, frame, msg, **kwargs):
        """
        Show a popup window with message
        :param frame: name of parent window
        :param msg: message to display
        :return:
            txt - text widget instance
            dlg - top level widget instance
        """
        if "height" in kwargs.keys():
            dlg = Toplevel()
        else:
            dlg = Toplevel()
        width=max([len(mm) for mm in msg.split('\n')])
        height=min([msg.count('\n')+1, 40])
        txt = Text(dlg, width=width, height=height)
        txt.insert("1.0", msg)
        txt.pack()
        if kwargs.has_key("location"):
            dx,dy=kwargs['location']
            w=txt.winfo_width()
            h=txt.winfo_height()
            x=frame.winfo_x()
            y=frame.winfo_y()
            loc="%dx%d+%d+%d" % (600,400, x+dx, y+dy)
            print "loc=",loc
            dlg.geometry(loc)
        if "wait" in kwargs.keys():
            dlg.focus_set()
            dlg.grab_set()
            dlg.transient(master=frame)
            dlg.wait_window(dlg)

        dlg.lift(aboveThis=frame)
        return txt, dlg

    def GetRootName(self, filename):
        basename = os.path.basename(filename)
        isuf = basename.find(".sbd")
        if isuf >= 0: basename = basename[0:isuf]
        return basename

    def SetupPreview(self):
        """Set range boxes, show preview of rt parameter variants, show help topic"""
        self.SetRangeBoxes()
        self.ShowPreview()
        self.HelpTopic()

    def SetRangeBoxes(self):
        """Set values of range boxes"""
        key = self.selectvariant.get()
        if hasattr(self, 'rangeskew'):
            self.rangeskew.delete(0, 'end')
            self.rangeskew.insert(0, 1.0)
        description,ranges = self.geninput.rtRange[key].split('$')
        self.description.delete(0,'end')
        self.description.insert(0,description.strip())
        self.rangestart.delete(0, 'end')
        self.rangefinis.delete(0, 'end')

        if ":" in ranges:
            r1, r2, skw = ranges.split(':')
            self.FixRangeStep=self.NiceStep(float(r2))
            self.rangestart.insert(0, r1)
            self.rangefinis.insert(0, r2)
            if hasattr(self, 'rangeskew'):
                self.rangeskew.delete(0,'end')
                self.rangeskew.insert(0, skw)
        else:
            self.rangestart.insert(0, ranges)
            self.rangefinis.insert(0, "")

    def ShowCoPreview(self, event):
        self.ShowPreview(covariant=1)

    def NiceStep(self, top):
        v=np.log10(abs(top))
        mantissa = int(v)
        fraction = v-mantissa
        step = 1
        if fraction < 0.30103:
            step = 0.1
        elif fraction < 0.69897:
            step = 0.2
        else:
            step = 0.5
        return step*10**mantissa

    def FixRange(self, event):
        """
        call back function for mousewheel events generated by range specification widgets
        :param event: a mousewheel event either up or down
        :return:
        """

        if event.widget == self.rangestart:
            #print "rangestart "
            # print self.rangestart.get()
            bot = float(self.rangestart.get())
            top = float(self.rangefinis.get())
            v=bot
            if event.delta < 0 and bot > 0:
                v-=self.FixRangeStep
            elif v < top-self.FixRangeStep:
                v+=self.FixRangeStep
            event.widget.delete(0,'end')
            event.widget.insert(0,"{:.3g}".format(v))

        if event.widget == self.rangefinis:
            bot = float(self.rangestart.get())
            top = float(self.rangefinis.get())
            v=top
            if event.delta > 0:
                v+=self.FixRangeStep
            elif v > bot+self.FixRangeStep:
                v-=self.FixRangeStep
            event.widget.delete(0,'end')
            event.widget.insert(0,"{:.3g}".format(v))

        if event.widget == self.rangesamples:
            v=int(self.rangesamples.get())
            if event.delta > 0 and v < 40:
                v+=1
            elif event.delta < 0 and v > 0:
                v-=1
            event.widget.delete(0,'end')
            event.widget.insert(0,"{}".format(v))

        if event.widget == self.rangeskew:
            v=float(self.rangeskew.get())
            if event.delta > 0 and v < 1:
                v+=0.1
            elif event.delta < 0 and v > 0:
                v-=0.1
            event.widget.delete(0,'end')
            event.widget.insert(0,"{:.1f}".format(v))
        self.ShowPreview()

    def ShowPreview(self, *args, **kwargs):
        key = self.selectvariant.get()
        skew = 1.0
        if hasattr(self, 'rangeskew'):
            skew = float(self.rangeskew.get())
        start = self.rangestart.get()
        finis = self.rangefinis.get()
        description,ranges = self.geninput.rtRange[key].split('$')
        #self.description.delete(0,'end')
        #self.description.insert(0,description)
        if ":" in ranges:
            ranges = "{}:{}".format(start, finis)  # in case user modifies start and finis from entry boxes
        number = int(self.rangesamples.get())
        lhs = self.SpreadValues(ranges, number, skew)
        cmd = "{}={}".format(key, lhs)
        if kwargs.has_key('covariant') and kwargs['covariant']:
            cmd += " &"
        cmd = self.geninput.DocString(cmd)
        self.PreviewLine(cmd, justify='left')

    def PreviewLine(self, cmd, **kwargs):
        """
        write stuff to preview line
        :param cmd:
        :param loc:
        :param kwargs: justify=-1, 0, 1  (left, center, right)
        :return:
        """
        self.previewbox.delete(0, 'end')
        self.previewbox.insert(0, cmd)
        if "justify" in kwargs.keys():
            self.previewbox.config(justify=kwargs['justify'])

    def SpreadValues(self, ranges, number, skew):
        """
        :param ranges: command string in format 'start:finish'
        :param number: integer number of values
        :param skew: skewness of distribution, small values are logarithmic, large values are linear
        :return: semicolon separated values
        """
        if ":" in ranges:
            r1, r2 = ranges.split(':')
            r1 = float(r1)
            r2 = float(r2)
            n = number
            if n == 1:
                lhs = r1
            elif skew == 1.0:
                values = np.linspace(r1, r2, n)
                vstr = ["{:.4g}".format(v) for v in values]
                lhs = ";".join(vstr)
            else:
                if r1 == 0 and skew == 0:
                    skew = 0.1
                dr = (r2 - r1) * skew * skew * skew
                values = np.logspace(np.log10(r1 + dr), np.log10(r2 + dr), n) - dr
                if r1 == 0:
                    values[0] = 0
                vstr = ["{:.4g}".format(v) for v in values]
                lhs = ";".join(vstr)
        elif ";" in ranges:
            nv = ranges.count(';') + 1
            iv = int(self.rangesamples.get()) % nv
            lhs = ranges.strip().split(';')[iv]

        return str(lhs)

    def AbortRun(self):
        self.abortit = True

    def NewRtParm(self, event):
        """Delete existing command listing"""
        self.sbdartoutput = ""
        self.caption.delete(1.0, 'end')
        self.diffbase = ''
        self.optionComparisonPlot.set(0)
        self.optionDiurnalPlot.set(0)
        self.PreviewLine("")
        #self.AddRtParm()

    def AddRtParm(self):
        parm = self.selectvariant.get()
        cmds = self.caption.get(1.0, 'end')
        newcmd = self.previewbox.get()
        values = newcmd.split('=')[1]

        if parm in newcmd:  # this test prevents copying bogus info messages
            if not cmds.startswith("#"):
                self.caption.insert("end", "# Description: \n#\n")

            if self.ParmNotFound(cmds, parm):
                self.caption.insert("end", newcmd + "\n")
            else:
                self.caption.delete(1.0,'end')
                self.caption.insert(1.0, self.ReplaceParm(cmds, newcmd))

            self.master.title("+ "+self.runname)
            root.update()
            self.runbtn.config(state="normal")
        else:
            self.PreviewLine("Select a parameter", justify='center')
            root.bell()

    def ParmNotFound(self, cmds, parm):
        lines = cmds.split('\n')
        for line in lines:
            if line.startswith('#'): continue
            if parm + "=" in line: return False
        return True

    def ReplaceParm(self, cmds, newcmd):
        rmcmd=newcmd.endswith('=') # no values indicate parameter should be removed
        parm=newcmd.split('=')[0]
        lines = cmds.split('\n')
        newcmds=""
        for line in lines:
            if line == '':
                continue
            elif line.startswith(parm + "="):
                if not rmcmd:
                    newcmds+=newcmd + "\n"
            else:
                newcmds+=line+"\n"
        return newcmds

    def SetCaption(self):
        self.caption.delete(1.0, 'end')
        fh = open("RUNS/" + self.menugroups[self.menugroup.get()], 'r')
        lines = fh.readlines()
        fh.close()
        for l in lines:
            self.caption.insert("end", l)

    def ClearPlot(self):
        if hasattr(self, 'ax'):
            self.ax.clear()

    def LoadFile(self, filename=""):
        """
        Open and read command file, set runname
        :return:
        """
        if not filename:
            filename = tkFileDialog.askopenfilename()
#                        filetypes=[('pkl files', '*.pkl'), ('sbd files', '*.sbd'), ('all files','*.*')])
        if filename:
            self.sbdartoutput = ""
            self.caption.delete('1.0', 'end')
            basename = self.GetRootName(filename)
            self.runbtn.config(state="disabled")
            try:
                if filename.endswith('sbd'):
                    self.load_flat_file(filename)
                else:
                    self.load_pickle_file(filename)
                self.runname = basename.split('.')[0]
                self.master.title(self.runname)
            except:
                self.caption.insert('1.0', "\n\n     File {} could not be read".format(filename))
        self.runbtn.config(state="normal")

    def load_flat_file(self, filename):
        '''
        load sbdart output in flat format
        :param basename:
        :param filename:
        :return:
        '''
        fh = open(filename, "r")
        sbdout = fh.read()
        fh.close()
        ind = sbdout.find("_DATA_")
        if ind < 0:
            self.caption.insert('1.0', sbdout)
        else:
            # self.addedparms = self.GetAddedParms(cmd[0: ind - 1])
            self.caption.insert('1.0', sbdout[0:ind - 1])
            root.update()
            self.sbdartoutput = sbdout[ind + 7:]
            if self.RtLoops(ingest=True):
                self.Plotit()
        if len(sbdout) < 300000:  # if its a short file, save filename for later recovery
            with open("RecentFile", 'w') as myfile:
                myfile.write(filename)

    def load_pickle_file(self, filename):
        '''
        load sbdart output in pickle format
        :param basename:
        :param filename:
        :return:
        '''
        with open(filename, 'rb') as fh:
            sbdout = pickle.load(fh)
            if 'COMMAND' in sbdout.keys():
                cmd = sbdout['COMMAND']
            else:
                print("keys in picklefile: ")
                for k in sbdout.keys():
                    print k
                    root.quit()

            self.caption.insert('1.0', cmd)
            root.update()
            self.yvariable=OrderedDict(sbdout)
            self.yrange=[float('inf'), float('-inf')]
            self.optionComparisonPlot.set(0)
            self.optionDiurnalPlot.set(0)
            self.diffbase = ''
            self.PreviewLine('')
            self.sbdartoutput = ""
            self.xvariable, self.xlabel = self.geninput.CycleSetup(cmd)
            self.parser = RtReader.RtReader(cmd, BooleanVar)
            if self.parser.IOUT in [20,21]:
                self.parser.zen = self.yvariable['ZEN']
                self.parser.phi = self.yvariable['PHI']
            self.SetRtMenu()
            self.SetGroupMenu()
            self.Plotit()
        if sys.getsizeof(sbdout) < 300000:  # if its a short file, save filename for later recovery
            with open("RecentFile", 'w') as myfile:
                myfile.write(filename)

    def WriteFile(self):
        """
        Write current command file to a sbd file
        :return:
        """
        fh = tkFileDialog.asksaveasfile(mode='w', initialdir = "RUNS", defaultextension=".sbd",
                                        filetypes=(("sbd files", "*.sbd"), ("All files", "*.*")), title='Write File')
        if fh:
            txt = self.caption.get('1.0', 'end')
            fh.writelines(txt)
            if self.sbdartoutput:
                fh.writelines('_DATA_\n')
                fh.writelines(self.sbdartoutput)
            fh.close()
            self.runname = self.GetRootName(fh.name)

    def PickleSave(self, mode = 'saveas'):
        """
        Write current command file to a sbd file
        :return:
        """

        if mode == 'saveas':
            fh = tkFileDialog.asksaveasfile(mode='wb', initialdir = 'RUNS', defaultextension=".pkl",
                                            filetypes=(("pkl files", "*.pkl"), ("All files", "*.*")), title='Save File')
            self.runname = str(os.path.basename(fh.name))
            if '.' in self.runname:
                self.runname = self.runname.split('.')[0]
        else:
            name =  "RUNS{}{}{}".format(os.path.sep, self.runname,'.pkl')
            fh = open(name, 'wb')

        if fh:
            cmd = self.caption.get('1.0', 'end')
            iout = int(self.geninput.IOUTformat)
            if iout == 10:
                obj = OrderedDict(self.yvariable)
                obj['COMMAND'] = cmd
                obj[self.xlabel] = np.array(self.xvariable).astype(float)
            elif iout in [20,21]:
                obj = OrderedDict(self.yvariable)
                obj['COMMAND'] = cmd
                obj['ZEN'] = self.parser.zen
                obj['PHI'] = self.parser.phi
            else:
                if iout == 11:
                    ky = 'ZZ'
                elif iout in [1,2]:
                    ky = 'WL'
                obj = OrderedDict()
                xvec = []
                for key in self.yvariable:
                    if key.startswith(ky):
                        if not xvec:
                            xvec = self.yvariable[key]
                    else:
                        obj[key]=self.yvariable[key]
                obj['COMMAND'] = cmd
                obj[ky]=xvec

            pickle.dump(obj, fh, -1)
            fh.close()
            self.master.title(self.runname)
            root.update()

    def ValidateCommands(self):
        """
        scan cmd text,
        make sure constant parms only appear after variable parms
        remove any leading or trailing white space
        add certain required constants parameters
        :return:
        """
        variableParms, constantParms, parms = self.GetParmsInCmd() # For clarity, add default switch values

        # add commands required by context

        if 'TBAER' in parms and not 'IAER' in parms:
            constantParms.append(self.geninput.DocString("IAER=1"))
        if 'TAERST' in parms:
            if not 'JAER' in parms:
                constantParms.append(self.geninput.DocString("JAER=1"))
            if not 'ZAER' in parms:
                constantParms.append("ZAER=15")
        if not 'IOUT' in parms:
            constantParms.append(self.geninput.DocString("IOUT=10"))

        self.caption.delete('1.0','end')
        rtcmds = ''
        if variableParms:
            rtcmds += '\n'.join(variableParms)+'\n'
        if constantParms:
            rtcmds += '\n'.join(constantParms)+'\n'
        self.caption.insert('1.0', rtcmds)

        # check for required input files

        requiredfiles = { 'IDATM=0':  'atms.dat',
                          'ISALB=-1': 'albedo.dat',
                          'NRE=0':    'usrcld.dat',
                          'IAER=-1':  'aerosol.dat',
                          'KDIST=-1': 'cktau.dat',
                          'ISAT=-1':  'filter.dat',
                          'NF=-1':    'solar.dat',
                          'NF=-2':    'cktau.dat'}
        msg=''


        for rf,fn in requiredfiles.iteritems():
            if (rf in constantParms) and not os.path.isfile(fn):
                msg+= "{} requires input file {}\n".format(rf,fn)

        nre='NRE=0'
        for c in variableParms:
            if c.startswith(nre):
                msg+= "{} requires input file {}\n".format(nre, requiredfiles(nre))

        if msg:
            self.Popup(self.framegraph, "Files not found:\n\n"+msg, wait=True)
            return False
        else:
            return True

    def GetParmsInCmd(self, **kwargs):
        '''
        scans command box text for variableParms, constantParms and parms list
        if keyword clean is true all comments are striped
        :return:
        variableParms   - list of commands that specify variable parameters
        constantParms   - list of commands that specify constant parameters
        parms           - string that contains ordered list of rt parameters
        '''
        txt = str(self.caption.get('1.0', 'end'))
        constantParms = []
        variableParms = []
        parms = ''

        if 'clean' in kwargs and kwargs['clean']:
            clean=True
        else:
            clean=False

        for pp in txt.split('\n'):
            p = pp.strip()
            if p.startswith('#'):
                if clean:
                    continue
                else:
                    variableParms.append(p)
                continue
            elif '=' in p:
                if not clean:
                    p = p.split('#')[0]
                pp = p.split('=')[0]
                parms += pp.upper() + ';'
                if ';' in p:
                    variableParms.append(p)
                else:
                    constantParms.append(p)

        return variableParms, constantParms, parms

    def SaveLastTry(self):
        """
        save last cmd text to ~RunRT.sbd
        :return:
        """
        txt = str(self.caption.get('1.0', 'end'))
        dir = ""
        if os.path.isdir("RUNS"):
            dir = "RUNS{}".format(os.path.sep)
        fh = open("{}~RunRT.sbd".format(dir), 'w')
        if "_DATA_" in txt:
            l = txt.index("_DATA_")
            fh.writelines(txt[0:l])
        else:
            fh.writelines(txt)
        fh.close()

    def WriteInput(self, inp):
        lun = open('./INPUT', 'w')
        lun.write("&INPUT\n")
        lun.write(inp)
        lun.write("/\n")
        lun.close()

    def ExecuteCmds(self, i, ingest):
        """
        Run sbdart for a iteration i
        :param i: iteration number
        :param ingest: read from sbdart output file if true, otherwise run sbdart
        :return: sbdart output for iteration i, input labels for iteration i
        """
        inp, lbls = self.geninput.CycleInput(i)
        niter=self.geninput.Niter()
        out = ""
        if ingest:
            buf = self.sbdartoutput.split('\n')
            nlines = len(buf)
            if nlines == niter:
                out = buf[i]
            else:
                nl = nlines/niter
                out='\n'.join(buf[i*nl:i*nl+nl])
        else:
            try:
                self.WriteInput(inp)
                proc = subprocess.Popen([self.sbdartexe], stdout=subprocess.PIPE)
                (out, err) = proc.communicate()
                self.sbdartoutput += out

            except:
                self.Popup(self.framegraph, err)
                print "out:", out
                print "err:", err

        return out, lbls

    def RtLoops(self, **kwargs):
        """
        run sbdart over command loops. changes xvariable and yvariable
        :return: True if all rt loops complete successfully, False otherwise
        """
        self.yrange=[float('inf'), float('-inf')]
        self.optionComparisonPlot.set(0)
        self.optionDiurnalPlot.set(0)
        self.diffbase = ''
        self.PreviewLine('')

        if kwargs.has_key('ingest') and kwargs['ingest']:
            ingest = True
        else:
            ingest = False
            self.sbdartoutput = ""
        cmds = self.caption.get("1.0", 'end')
        if cmds.count('=') == 0:
            return False
        self.xvariable, self.xlabel = self.geninput.CycleSetup(cmds)
        if 'Error' in self.xlabel:
            self.Popup(self.framegraph, self.xlabel)
            return False
        self.parser = RtReader.RtReader(cmds, BooleanVar)
        self.yvariable = OrderedDict()
        labels = []
        niter = self.geninput.Niter()
        reporttime = 1
        timemark = timeit.default_timer() + reporttime
        for i in range(0, niter):
            out, lbls = self.ExecuteCmds(i, ingest)
            if self.parser.IOUT == 10:
                if len(lbls) > 0:
                    label = " ".join(sorted(lbls[1:]))
                else:
                    label=""
            elif len(lbls) > 0:
                label = " ".join(sorted(lbls))
            else:
                label = ""
            labels.append(label)
            fraction = "{:.0f}%".format(i * 100.0 / niter)
            if ingest:
                timenow = timeit.default_timer()
                if timenow > timemark:
                    self.PreviewLine(fraction, justify='center')
                    root.update()
                    timemark = timenow + reporttime
            else:
                wpv=100
                updatelabel = label + " "*(wpv-len(label)-30)+fraction+" "*30
                self.PreviewLine(updatelabel, justify='right')
                root.update()

            self.parser.RT(out, label, self.yvariable)
            if self.abortit:
                return False

        self.PreviewLine("")
        self.SetRtMenu()
        self.SetGroupMenu()
        return True



    def RunSBDART(self):
        for fn in glob.glob("SBDART_WARNING*"):
            os.remove(fn)
        self.abortit=False
        self.abortbtn.pack()
        self.addvariant.pack_forget()
        self.runbtn.config(state="disabled")
        if self.ValidateCommands():

            self.SaveLastTry()
            try:
                if self.RtLoops():
                    self.Plotit()
                    self.master.title("* "+self.runname)
                    root.update()

            except:
                self.PreviewLine("No plots available.", justify='center')
                self.ViewOutput()

        self.abortbtn.pack_forget()
        self.addvariant.pack()
        self.runbtn.config(state="normal")
        for fn in glob.glob("SBDART_WARNING*"):
            self.ShowFilePopup(fn)



if len(sys.argv) == 1:
    matplotlib.use('Tkagg')
    root = Tk()
    my_gui = RunRT(root)
    root.mainloop()
else:                           # if command line supplies filenames
    newdir = ""
    for fn in sys.argv[1:]:
        if fn.endswith(os.sep):
            newdir=fn
            continue
        with open(fn, 'r') as myfile:
            cmds=myfile.read()
        ind = cmds.find("_DATA_")
        sbdartoutput=""
        cmds=cmds[0:ind]
        gi=GenInput.GenInput()
        junk, msg = gi.CycleSetup(cmds)
        if 'Error' in msg:
            print msg
            exit(1)

        for i in range(0, gi.Niter()):
            inp, lbls = gi.CycleInput(i)
            with open("INPUT",'w') as inpfile:
                inpfile.write("&INPUT\n{}\n/\n".format(inp))
            'fn {}'.format(fn)
            proc = subprocess.Popen([sbdartexe], stdout=subprocess.PIPE)
            (out, err) = proc.communicate()
            sbdartoutput += out
        if newdir:
            if not os.path.exists(newdir):
                os.makedirs(newdir)
            sfn = newdir + os.path.basename(fn)
            os.rename(fn, sfn)
        with open(fn, 'w') as myfile:
            myfile.write("{}\n_DATA_\n{}".format(cmds,sbdartoutput))
