"""
Various functions and classes to simplify dynamics displays for class
"""

import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import time


                
def right_hand_fig(k, extra_width=0, clearfig=True):
    """
    Makes a figure on the right of a Mac screen.
    
    PARAMS:
        
        k   number of figure; if 1, fig will be placed at the top, if 2, at the bottom
        
        extra_width=0  Default width is 640; this adds to it and extends the figure leftwards.
        
        clearfig=True  whether to clear the figure after creating it
        
    RETURNS:
    
        fig   The created figure object
    """
    fig = plt.figure(k)
    if clearfig:
        plt.clf()
        
    manager = plt.get_current_fig_manager()
    wpos = list(manager.window.geometry().getRect())
    wpos[0] = 1270 - extra_width; 
    wpos[2] = 640 + extra_width
    wpos[3] = 550
    if k==1:
        wpos[1] = 50
    elif k==2:
        wpos[1] = 625
        
    manager.window.setGeometry(wpos[0], wpos[1], wpos[2], wpos[3])
    manager.window.raise_()
    
    fig.canvas.draw()
    fig.canvas.flush_events()
    return fig


class kb_monitor:    
    """
    An instance of this class sets up monitoring of keystrokes and button presses in a figure.
    
    Create as kb_monitor(fig) where fig is a matplotlib.pyplot figure handle; after this, calling 
    the methods keylist() or buttonlist() will return lists of the keys and the buttons pressed, respectively.
    """
    def __init__(self, fig):
        """
        Create and return an instance of a key and button monitor, attached to a specific figure
        
        PARAMS:
            fig   a handle to a matplotlib figure, that this object will monitor
            
        RETURNS:
            none
        """
        self.__myfig = fig
        self.__keys_pressed = []
        self.__buttons_pressed = []
        self.__bcid = fig.canvas.mpl_connect('button_press_event', self.__button_callback)        
        self.__kcid = fig.canvas.mpl_connect('key_press_event', self.__key_callback)        

    def __button_callback(self, event):
        self.__buttons_pressed += [(event.inaxes, event.xdata, event.ydata)]

    def __key_callback(self, event):
        self.__keys_pressed += [event.key]

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
        
    