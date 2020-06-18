try:
    from tkinter import Tk, Toplevel, Text, Entry, Button, Frame, PhotoImage
except:
    from Tkinter import Tk, Toplevel, Text, Entry, Button, Frame, PhotoImage
import sys
from pathlib import Path
import os


class Docview:

    def __init__(self, master, fn, **kwargs):

        self.file = fn
        self.widget = None
        self.text = Text()
        self.frame = master
        self.location = kwargs.get('location', None)
        self.wait = kwargs.get('wait', None)
        self.find = kwargs.get('find', None)

        if getattr(sys, 'frozen', False):
            folder = Path(sys._MEIPASS)
        else:
            folder = Path(__file__).parent
        self.delete_img = PhotoImage(file=f'{folder}{os.path.sep}delete.pbm')
        self.up_img = PhotoImage(file=f'{folder}{os.path.sep}uparrow.pbm')
        self.down_img = PhotoImage(file=f'{folder}{os.path.sep}downarrow.pbm')
        self.upbutn = Button()
        self.dnbutn = Button()
        self.deletebutn = Button()
        self.entry = None
        self.doc = ''


    def Search(self, pattern):
        '''
        normal top-down search
        :param pattern:
        :return:
        '''
        try:
            self.text.tag_remove('hilight','1.0','end')
            ind = self.text.search(pattern, '1.0', stopindex="end", nocase=True)
            self.text.mark_set("insert", ind)
            self.text.see("insert")
            nchar = len(pattern)
            self.text.tag_add('hilight', ind,  ind+"+{}c".format(nchar))

        except:
            pass

    def RecursiveSearch(self, pattern, upward=False):
        '''
        Search forward from current position.  Search backward if upward is True.
        :param pattern: search pattern
        :param upward: reverse search if true
        :return: Null
        '''
        if len(pattern)==0:
            ind = '1.0' if upward else 'end'
            self.text.mark_set("insert", ind)
            self.text.see("insert")
            return
        try:

            insert = self.text.index("insert")
            if upward:
                point = insert + "-1c"
                ind = self.text.search(pattern, point, stopindex="1.0", backwards=True, nocase=True)
            else:
                point = insert + "+1c"
                ind = self.text.search(pattern, point, stopindex="end", nocase=True)

            self.text.mark_set("insert", ind)
            self.text.see("insert")
            self.text.tag_remove('hilight','1.0','end')
            nchar = len(pattern)
            self.text.tag_add('hilight', ind,  ind+"+{}c".format(nchar))
            self.text.update()
        except:
            pass

    def SearchUp(self):
        pattern = str(self.entry.get())
        self.RecursiveSearch(pattern, True)

    def SearchDn(self):
        pattern = str(self.entry.get())
        self.RecursiveSearch(pattern, False)

    def FindEntry(self, event):
        '''triggered by keyrelease event in find entry box'''
        pattern = str(self.entry.get())
        nchar = len(pattern)
        self.text.tag_remove('hilight','1.0','end')
        if nchar == 0:
            return
        try:
            ind = self.text.search(pattern, "1.0", stopindex="end", nocase=True)
            self.text.mark_set("insert", ind)
            self.text.see("insert")
            self.text.tag_add('hilight', ind,  ind+"+{}c".format(nchar))
            self.text.update()
        except:
            pass

    def DeleteSearch(self):
        self.entry.delete(0,'end')
        self.text.tag_remove('hilight', '1.0','end')

    def Popup(self):

        if not self.doc:
            with open(self.file, 'r') as fh:
                self.doc = fh.read()

        if self.widget:     # widget already open
            if self.frame:
                self.widget.lift(aboveThis=self.frame)
            return

        #print 'length of doc ', len(self.doc)

        self.widget = Toplevel()
        if self.find:
            self.sframe = Frame(self.widget)
            self.sframe.pack()
            self.upbutn = Button(self.sframe,width=17, command = self.SearchUp, image=self.up_img)
            self.upbutn.pack(side='left')
            self.dnbutn = Button(self.sframe,width=17, command = self.SearchDn, image=self.down_img)
            self.dnbutn.pack(side='left')
            self.entry = Entry(self.sframe, width=50)
            self.entry.bind('<KeyRelease>', self.FindEntry)
            self.entry.pack(side='left')
            self.deletebutn = Button(self.sframe, width=17, command = self.DeleteSearch)
            self.deletebutn.config(image=self.delete_img)
            self.deletebutn.pack(side='left')

        width=max([len(line) for line in self.doc.split('\n')])
        height=min([self.doc.count('\n')+1, 40])
        self.text = Text(self.widget, width=width, height=height)
        self.text.tag = self.text.tag_configure('hilight', background='#ffff00')
        self.text.insert("1.0", self.doc)
        self.text.pack()

        if self.location and self.frame:
            dx,dy=self.location
            w=self.text.winfo_width()
            h=self.text.winfo_height()
            x=self.frame.winfo_x()
            y=self.frame.winfo_y()
            loc="%dx%d+%d+%d" % (600,400, x+dx, y+dy)
            self.widget.geometry(loc)

        if self.wait:
            self.widget.focus_set()
            self.widget.grab_set()
            if self.frame:
                self.widget.transient(master=self.frame)
            self.widget.wait_window(self.widget)

        if self.frame:
            self.widget.lift(aboveThis=self.frame)

        self.widget.protocol("WM_DELETE_WINDOW", self.Kill)

    def Kill(self):
        self.widget.destroy()
        self.widget = None

if __name__ == "__main__":
    root = Tk()
    dv = Docview(root, 'rtdoc.txt', find=True)
    dv.Popup()
    root.mainloop()
