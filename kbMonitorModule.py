"""
kbMonitorModule.py -- This module contains a single class, kb_monitor.

To get documentation, ask for help on kbMonitorModule.kb_monitor.

"""

        

class kb_monitor:    
    """
    An instance of this class sets up monitoring of keystrokes and button presses in a figure.
    
    Create as kb_monitor(fig) where fig is a matplotlib.pyplot figure handle; after this, calling 
    the methods keylist() or buttonlist() will return lists of the keys and the buttons pressed,
    respectively.
    
    You can also create it as kb_monitor(fig, callback=func), in which case func() is called
    immediately after any button or key press; within func() you can then call keylist() or
    buttonlist() to see events. func() should take one obligatory argument, which will be
    the kb_monitor object.
    
    From Julia, you would furst use PyCall to load the module
    
        using PyCall
        # The following line is PyCall-ese for "add the current directory to the Python path"
        unshift!(PyVector(pyimport("sys")["path"]), "")
        # We use Python to enable callbacks from the figures.
        @pyimport kbMonitorModule
        
    And then you would use it with the PyCall syntax, for example
    
        function my_callback_func(BP)
            print(BP[:buttonlist]()); print("\n")
        end
        
        using PyPlot
        pygui(true)
        BP = kbMonitorModule.kb_monitor(figure(1), callback=my_callback_func)

    """
    def __init__(self, fig, callback=None):
        """
        Create and return an instance of a key and button monitor, attached to a specific figure
        
        PARAMS:
            fig   a handle to a matplotlib figure, that this object will monitor
            
        RETURNS:
            none
        """
        self.__myfig = fig
        self.__user_callback = callback
        self.__keys_pressed = []
        self.__buttons_pressed = []
        self.__bcid = fig.canvas.mpl_connect('button_press_event', self.__button_callback)        
        self.__kcid = fig.canvas.mpl_connect('key_press_event', self.__key_callback) 

    def __button_callback(self, event):
        self.__buttons_pressed += [(event.inaxes, event.xdata, event.ydata)]
        if self.__user_callback != None:
            self.__user_callback(self)

            
    def __key_callback(self, event):
        self.__keys_pressed += [event.key]
        if self.__user_callback != None:
            self.__user_callback(self)


    def figure_handle(self):
        """
            PARAMS:
                None
            RETURNS:
                the handle to the figure associated with the kb_monitor               
        """
        return self.__myfig
    
        
    def keylist(self):
        """
            PARAMS:
                None
            RETURNS:
                A list of single characters, indicating keys pressed, in order of their pressing.
                Empty list means no keys pressed so far.
        """
        return self.__keys_pressed
    
    def buttonlist(self):
        """
            PARAMS:
                None
            RETURNS:
                A list of tuples, each containing three elements. The first element is an axis handle
                (if the button was pressed while the mouse was pointing within an axis, otherwise None); 
                the next two are the x position and y position of the mouse, in data coordinates (None and
                None is the mouse was not pointing within an axis)
        """
        return self.__buttons_pressed
        
    def clear_keylist(self):
        self.__keys_pressed = []
        
    def clear_buttonlist(self):
        self.__buttons_pressed = []
        
    def __del__(self):
        self.__myfig.canvas.mpl_disconnect(self.__bcid)
        self.__myfig.canvas.mpl_disconnect(self.__kcid)
        
    