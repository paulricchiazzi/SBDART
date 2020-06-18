try:
    from tkinter import Entry
except:
    from Tkinter import Entry

class Spinner(Entry):
    """
    Customized Entry widget to allow simple selection from a list of entries
    By default, the return key is bound to a search function, but this can be over-ridden
    keywords:
        Return  -- command that over-rides search function
        values  -- a list of values
        init    -- initial value in values
        command -- a command to run each time selection changes
        wrap    -- if true, roll up or down events wrap around at top and bottom of list
    """

    def __init__(self, master, **kwargs):

        self.values=[]
        edict = {}
        init = ""
        for k in kwargs.keys():
            if not k in ('Return', 'values', 'init', 'command','wrap'):
                edict[k]=kwargs[k]
        # edict['validate']='key'
        # okaycmd=self.register(self.isOkay)
        # edict['validatecommand']=(okaycmd, '%S')


        Entry.__init__(self, master, **edict)


        if 'Return' in kwargs:
            self.bind('<Return>', kwargs['Return'])
        else:
            self.bind('<Return>', self.Search)

        self.wrap=False

        if 'values' in kwargs:
            self.values = kwargs['values']
        if 'init' in kwargs:
            init=kwargs['init']
        if 'command' in kwargs:
            self.CommandCallBack=kwargs['command']
        if 'wrap' in kwargs:
            self.wrap=True

        if self.values:
            self.Index=0
            if init and init in self.values:
                self.Index=self.values.index(init)
            self.insert(0, self.values[self.Index])

        self.bind('<MouseWheel>', self.MouseWheelEvent)
        self.bind('<Up>', self.FlipUp)
        self.bind('<Down>', self.FlipDn)

    def MouseWheelEvent(self, event):
        self.FlipToValue(event.delta)

    def FlipUp(self, event):
        event.char=None
        self.FlipToValue(2)

    def FlipDn(self, event):
        event.char=None
        self.FlipToValue(-2)

    def isOkay(self, S):
        if ord(S)>= ord(' ') and ord(S) <= ord('~'):
            return True
        else:
            return False


    def FlipToValue(self, delta):
        n=len(self.values)
        if n == 0:
            return
        if delta > 0 and (self.Index < n - 1 or self.wrap):
            self.Index = (self.Index + 1) % n
        elif delta < 0 and (self.Index > 0 or self.wrap):
            self.Index = (self.Index + n - 1) % n
        self.delete(0, 'end')
        self.insert(0, self.values[self.Index])
        if self.CommandCallBack:
            self.CommandCallBack()

    def Search(self, event):
        name = str(self.get()).upper()
        for nm in self.values:
            rtnm = str(nm).upper()
            if rtnm.startswith(name):
                self.Index = self.values.index(nm)
                self.delete(0,'end')
                self.insert(0,self.values[self.Index])
                if self.CommandCallBack:
                    self.CommandCallBack()
