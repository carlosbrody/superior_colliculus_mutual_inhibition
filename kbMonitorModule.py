"""
kbMonitorModule.py -- This module contains a single class, kb_monitor.

To get documentation, ask for help on kbMonitorModule.kb_monitor.

To use in Julia, do
    using PyCall
    if PyVector(pyimport("sys")["path"])[1] != ""
        # Then the following line is PyCall-ese for "add the current directory to the Python path"
        unshift!(PyVector(pyimport("sys")["path"]), "")
    end

    @pyimport kbMonitorModule
    
functions will then be available as kbMonitorModule.funcname()

"""

import matplotlib.widgets as Wid


def text_box(ax, leftside_label, initial_text, user_callback=None):
    """
    tbox = text_box(ax, leftside_label, initial_text, user_callback=None)

Puts up, in the passed axes ax, a text box with leftside_label (a string)
to its left. Initial text in the box will initial_text (also a string).
If the user_callback is set,, that function will be called with one 
parameter, the string in the box.

    You can also access the string currently in the box through
    
    tbox.text
        
Or, in Julia, with the PyCall syntax

    rab[:text]
    
And you can set the value with 

    rab[:set_value] = your_string

    """
    tbox = Wid.TextBox(ax, leftside_label, initial=initial_text)
    if user_callback != None:
        tbox.on_submit(user_callback)

    return tbox



def check_buttons(rax, labels, actives, user_callback=None):
    """
    rab = check_buttons(ax, labels, actives, user_callback=None)

Puts up check buttons (any number one selected at a time) for the strings
passed in the vector of strings "labels", in the matplotlib axes given by "ax".  
actives must be a list of Booleans, same length as labels,
indicating which buttons should start as ON.
On click, if the function user_callback is set, it will be called, with a string 
as its single parameter-- that string will be the label of the clicked button. 

You can also access which buttons are currently on through
    
    rab.get_status()
        
Or, in Julia, with the PyCall syntax

    rab[:get_status]()
    
which returns a tuple of Booleans, with true for those that are ON,
false for those that aren't.  A list of string labels can also be obtained 
through

    rab[:labels]
    
    """
    checks = Wid.CheckButtons(rax, labels, actives)
    if user_callback != None:
        checks.on_clicked(user_callback)

    return checks
    

def radio_buttons(rax, labels, user_callback=None):
    """
    rab = radio_buttons(rax, labels, user_callback=None)

Puts up radio buttons (only one selected at a time) for the strings
passed in the vector of strings "labels", in the matplotlib axes given by "rax".  
On click, if the function user_callback is set, it will be called, with a string 
as its single parameter-- that string will be the label of the selected button. 

You can also access the string of the currently selected button through
    
    rab.value_selected
        
Or, in Julia, with the PyCall syntax

    rab[:value_selected]
     
    """
    radio = Wid.RadioButtons(rax, labels)
    if user_callback != None:
        radio.on_clicked(user_callback)

    return radio
    

    
class kb_monitor:    
    """
    BP = kbMonitorModule.kb_monitor(fig, callback=None, userData=None)
    
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
    def __init__(self, fig, callback=None, userData=None):
        """
        Create and return an instance of a key and button monitor, attached to a specific figure
        
        PARAMS:
            fig   a handle to a matplotlib figure, that this object will monitor
            
        RETURNS:
            none
        """
        self.__myfig = fig
        self.__user_callback = callback
        self.__user_data     = userData
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
        """
            PARAMS:
                None
            RETURNS:
                None -- but clears the internal list of keys pressed               
        """        
        self.__keys_pressed = []
        
    def clear_buttonlist(self):
        """
            PARAMS:
                None
            RETURNS:
                None -- but clears the internal list of buttons clicked               
        """        
        self.__buttons_pressed = []
        
    def set_userdata(self, userData):
        """
            PARAMS:
                userData, information to be stored internally
            RETURNS:
                None
        """
        self.__user_data = userData

    def get_userdata(self):
        """
            PARAMS:
                None
            RETURNS:
                whatever was stored internally as the user Data.
        """
        return self.__user_data
        
        
    def __del__(self):
        self.__myfig.canvas.mpl_disconnect(self.__bcid)
        self.__myfig.canvas.mpl_disconnect(self.__kcid)
        
    