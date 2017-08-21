from Tkinter import Tk, Toplevel, Text, Entry, Label

class Docview:

    def __init__(self, master, fn, **kwargs):

        self.file = fn
        self.widget = None
        self.text = Text()
        self.visible = False
        self.frame = master
        self.location = kwargs.get('location', None)
        self.wait = kwargs.get('wait', None)
        self.doc = ''

    def Search(self, pattern):
        ind = self.text.search(pattern, "1.0", stopindex="end")
        self.text.mark_set("insert", ind)
        self.text.see("insert")

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

        width=max([len(line) for line in self.doc.split('\n')])
        height=min([self.doc.count('\n')+1, 40])
        self.text = Text(self.widget, width=width, height=height)
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
    dv = Docview(root, 'rtdoc.txt')
    dv.Popup()
    root.mainloop()
