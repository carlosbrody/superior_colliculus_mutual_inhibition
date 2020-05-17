module GeneralUtils

export replacer

using Printf

"""
    tbin(tvector, t)

Returns indmin(abs.(tvector-t)) -- just a shorthand way
for finding, in a vector of time bins tvector, the bin
that corresponds the closest to time t.
"""
function tbin(tvector, t)
    return indmin(abs.(tvector-t))
end


"""
    newV = vstack_and_NaN_pad(V; ntrials=1)

if V is a vector, then makes a column vector that is ntrials copies
of V, stacked on top of each other; the copies are separated by one
element, containing NaN.

if V is m-by-n-by-k, then does not stack copies. Instead it takes
each V[:,:,i], takes its transpose, and then stacks all of those
vertically, again with a layer of NaN between each stack.

The net result of this is that if you have a time vector t
and some data V that is nunits x length(t) x ntrials long,
you can do

    `plot(vstack_and_NaN_pad(t, size(V,3)), vstack_and_NaN_pad(V))`

and all the trials will get plotted all at once, with each unit being a
unique color. nunits line handles will get returned by `plot()`.
"""
function vstack_and_NaN_pad(oV; ntrials=1)
    if length(size(oV))==3
        ntrials = size(oV,3)
        lent    = size(oV,2)
        V = Array{Float64}(lent*ntrials + ntrials-1, size(oV,1))
        for i=1:ntrials
            offset = (i-1)*(lent+1)
            V[(1:lent)+offset,:] = oV[:,:,i]'

            if i<ntrials
                V[lent+offset+1,:] = NaN
            end
        end
    else
        # @printf("size(oV)=\n"); print(size(oV)); print("\n")
        lent    = length(oV)
        # @printf("lent=%d ntrials=%d lent*ntrials + ntrials-1=%d\n", lent, ntrials, lent*ntrials + ntrials-1)
        V = Array{Float64}(lent*ntrials + ntrials-1)
        for i=1:ntrials
            offset = (i-1)*(lent+1)
            V[(1:lent)+offset] = oV[:]'

            if i<ntrials
                V[lent+offset+1] = NaN
            end
        end
    end
    return V
end


# DON'T MODIFY THIS FILE -- the source is in file General Utilities.ipynb. Look there for further documentation and examples of running the code.






"""
    append_to_file(filename, str)

Opens filename, appends str to it, and closes filename.

If filename is not a string but us type Base.PipeEndpoint (an IO stream)
then simply prints to it, without trying to open or close
"""
function append_to_file(filename, str)
    if typeof(filename)<:IO
        @printf(filename, "%s", str)
    else
        fstr = open(filename, "a")
        @printf(fstr, "%s", str)
        close(fstr)
    end
end



"""
evaluated_expression = replacer(P::String, mypars)

Given a string representing an expression to be evaluated, and a dictionary of
keys that are strings representing variables with values that are numbers,
parses the string into an expression, and substitutes any variables matching
the keys in mypars with the corresponding numeric values; finally, evaluates the
expression and returns the result.

# PARAMETERS:

- P::String   The expression to be evaluated, for example "t1*10 + sqrt(t2)"

- mypars      A dictionary mapping variable names to values, for example
                Dict("t1"=>5, "t2"=>100)

# RETURNS:

The result of evaluating the corresponding expression. Any variables that cannot
be instantiated into values will result in an Undefvar error

# EXAMPLE:

```jldoctest
julia>  replacer("t1*10+sqrt(t2)", Dict("t"=>-3, "t1"=>5, "t2"=>100))

15

```

"""
function replacer(P::String, mypars)
    return replacer(Meta.parse(P), mypars)
end


"""
evaluated_expression = replacer(P::Expr, mypars)

Given an expression to be evaluated, and a dictionary of
keys that are strings representing variables with values that are numbers,
substitutes any variables matching
the keys in mypars with the corresponding numeric values; finally, evaluates the
expression and returns the result.

# PARAMETERS:

- P::String   The expression to be evaluated, for example parse("t1*10 + sqrt(t2))"

- mypars      A dictionary mapping variable names to values, for example
                Dict("t1"=>5, "t2"=>100)

# RETURNS:

The result of evaluating the corresponding expression. Any variables that cannot
be instantiated into values will result in an Undefvar error

# EXAMPLE:

```jldoctest
julia>  replacer(parse("t1*10+sqrt(t2))", Dict("t"=>-3, "t1"=>5, "t2"=>100))

15

```

"""
function replacer(P::Expr, mypars)   # run through an expression tree, replacing known symbols with their values, then evaluate
    mypars = symbol_key_ize(mypars)
    ks = collect(keys(mypars))

    # otherwise, see if there are subarguments that we should work on
    if any(fieldnames(P).==:args)
        for i=1:length(P.args)
            # Not sure why, but need this first check for a tuple Expr to deal with fun.() syntax
            if typeof(P.args[i])<:Expr && P.args[i].head == :tuple
                for j=1:length(P.args[i].args)
                    P.args[i].args[j] = replacer(P.args[i].args[j], mypars)
                end
            elseif typeof(P.args[i])<:Symbol || typeof(P.args[i])<:Expr # && P.args[i].head != :tuple)
                # if you have a Symbol or Expr, go recursively into it
                P.args[i] = replacer(P.args[i], mypars)
            end
        end
    end
    # @printf("\nP = \n"); dump(P)
    # @printf("\nreplacing with \n"); dump(eval(P))
    return eval(P)
end


function replacer(P::Symbol, mypars)
    # if its a Symbol, see if it is in our dictionary, in which case replace it with its value
    idx = find(ks .== P)
    if length(idx)>0
        P = mypars[ks[idx[1]]]
    end
    return P
end

function replacer(P, mypars)
    return P
end

# DON'T MODIFY THIS FILE -- the source is in file General Utilities.ipynb. Look there for further documentation and examples of running the code.



"""
function print_vector(vec)

Takes a vector and uses @printf to put it on the screen with [%.3f, %.3f] format.

If passed a symbol (which must evaluate to a vector), then prints the string for that symbol,
an equals sign, the vector, and ends by adding a carriage return \n.
"""
function print_vector(vec)
    print_vector(STDOUT, vec)
end

"""
function print_vector(fname::String, vec)

Takes a vector and uses @printf to append it to file fname with [%.3f, %.3f] format.

If passed a symbol (which must evaluate to a vector), then prints the string for that symbol,
an equals sign, the vector, and ends by adding a carriage return \n.
"""
function print_vector(fname::String, vec)
    ostream = open(fname, "a")
    print_vector(ostream, vec)
    close(ostream)
end


"""
function print_vector(stream::IO, vec)

Takes a vector and uses @sprintf to put it on stream IO with [%.3f, %.3f] format.

If passed a symbol (which must evaluate to a vector), then prints the string for that symbol,
an equals sign, the vector, and ends by adding a carriage return \n.
"""
function print_vector(stream::IO, vec)

    if typeof(vec)==Symbol
        mystr = string(vec)
        @printf(stream, "%s = ", mystr);
        print_vector(stream, eval(vec))
        @printf(stream, "\n");
        return
    end

    @printf stream "["
    for p in [1:length(vec);]
        @printf(stream, "%.3f", vec[p])
        if p < length(vec) @printf(stream, ", "); end
    end
    @printf(stream, "]")
end


"""
function print_vector_g(vec)

Takes a vector and uses @printf to put it on the screen with [%g, %g] format.

If passed a symbol (which must evaluate to a vector), then prints the string for that symbol,
an equals sign, the vector, and ends by adding a carriage return \n.
"""
function print_vector_g(vec)
    print_vector_g(STDOUT, vec)
end


"""
function print_vector_g(fname::String, vec)

Takes a vector and uses @printf to append it to file fname with [%g, %g] format.

If passed a symbol (which must evaluate to a vector), then prints the string for that symbol,
an equals sign, the vector, and ends by adding a carriage return \n.
"""
function print_vector_g(fname::String, vec)
    ostream = open(fname, "a")
    print_vector_g(ostream, vec)
    close(ostream)
end


"""
function print_vector_g(stream::IO, vec)

Takes a vector and uses @printf to put it on stream with [%g, %g] format.

If passed a symbol (which must evaluate to a vector), then prints the string for that symbol,
an equals sign, the vector, and ends by adding a carriage return \n.
"""
function print_vector_g(stream::IO, vec)

    if typeof(vec)==Symbol
        mystr = string(vec)
        @printf(stream, "%s = ", mystr);
        print_vector_g(stream, eval(vec))
        @printf(stream, "\n");
        return
    end

    @printf stream "["
    for p in [1:length(vec);]
        @printf(stream, "%g", vec[p])
        if p < length(vec) @printf(stream, ", "); end
    end
    @printf(stream, "]")
end


"""
y = two_level_copy(x)

Like copy(x), but can go down a level. Can handle both Arrays and Dicts, otherwise gets confused.

EXAMPLE:

p = [1, 2, 3]
z = [p, 4]
y = Dict(:my=>p)

c = copy(z)
d = copy(y)

alpha = two_level_copy(z)
beta  = two_level_copy(y)
p[1]=1000

print("The inner levels of c and d are affected by the change to p:\n")
print(c); print("\n")
print(d); print("\n")
print("But the inner levels of alpha and beta are not:\n")
print(alpha); print("\n")
print(beta); print("\n")

"""
function two_level_copy(x)
    if typeof(x)<:Array
        y = copy(x)
        for i=1:length(x)
            if typeof(x[i])<:Tuple; y[i]=x[i];
            else y[i] = copy(x[i]) end;
        end
    elseif typeof(x)<:Dict
        y = copy(x)
        allkeys = keys(x)
        for i in allkeys
            if typeof(x[i])<:Tuple; y[i]=x[i];
            else y[i] = copy(x[i]) end;
        end
    else
        error(@sprintf("two_level_copy: Don't know how to handle type %s\n", typeof(x)))
    end
    return y
end




"""
new_fname = next_file(fbasename, ndigits)

Returns a numbered and presumably unused filename starting with the string fbasename, followed by an integer
digit. The returned integer will be one higher than the number of existing filenames starting with fbasename,
and will be written with ndigits numbers, using leading zeros if necessary.

# EXAMPLE:

If there are already 8 files starting with "Mydir/model" then

> next_file("Mydir/model_", 4)

"Mydir/model_0009"
"""
function next_file(fbasename, ndigits)
    mydir  = dirname(fbasename)
    myfile = basename(fbasename)
    if length(mydir)>0
        fnames = readdir(mydir)
    else
        fnames = readdir()
    end
    matched_filenames = Array{Bool}(length(fnames))
    for i=1:length(fnames)
        matched_filenames[i] = ismatch(Regex(@sprintf("^%s", myfile)), fnames[i])
    end

    mynum = length(find(matched_filenames))+1
    myname = @sprintf("%d", mynum)
    while length(myname)<ndigits
        myname = "0" * myname
    end

    if length(mydir)>0
        return mydir * "/" * myfile * myname
    else
        return myfile * myname
    end
end


"""
    fstring = num2fixed_string(n, ndigits)

Returns a string version of a positive integer, with
however many leading zeros are necessary to have
ndigits characters.

"""
function num2fixed_string(n, ndigits)
    if ndigits<=0
        error("ndigits must be bigger than zero")
    end

    if n<0
        error("n must be positive")
    end

    myname = @sprintf("%d", n)
    while length(myname)<ndigits
        myname = "0"*myname
    end

    return myname
end


# DON'T MODIFY THIS FILE -- the source is in file General Utilities.ipynb. Look there for further documentation and examples of running the code.



"""
ad = ascii_key_ize(d)

Given a dictionary that has keys that can be converted to strings, returns a copy with all
keys converted to strings
"""
function ascii_key_ize(d)
    ad = Dict()
    for k in keys(d)
        get!(ad, string(k), d[k])
    end
    return ad
end



"""
sd = symbol_key_ize(d)

Given a dictionary that has keys that can be converted to Symbols, returns a copy with all
keys converted to Symbols
"""
function symbol_key_ize(d)
    sd = Dict()
    for k in keys(d)
        get!(sd, Symbol(k), d[k])
    end
    return sd
end


"""
vks = vectorize_dict(dictionary, ks)

Given a dictionary (in which) all keys are either strings or Symbols, and all values are Float64s),
and an array ks of keys into that dictionary, returns a Float64 array the same size as ks containing
the values. Each key is checked as either itself or the string version of itself or the Symbol version
of itself.

Thus the following all return the same

a = Dict(:this=>33.4, "that"=>28.7)

vectorize_dict(a, ["this", "that"])
vectorize_dict(a, [:this, "that"])
vectorize_dict(a, [:this, :that])
"""
function vectorize_dict(dictionary, ks)
    output = Array{Float64}(size(ks))
    for i=1:length(ks)
        if haskey(dictionary, ks[i])
            output[i] = dictionary[ks[i]]
        elseif typeof(ks[i])<:Symbol && haskey(dictionary, string(ks[i]))
            output[i] = dictionary[string(ks[i])]
        elseif typeof(ks[i])<:String && haskey(dictionary, Symbol(ks[i]))
            output[i] = dictionary[Symbol(ks[i])]
        else
            print("Troublesome key: "); print(ks[i]); print("\n")
            error("Found neither key nor string(key) nor Symbol(key) in the dictionary")
        end
    end
    return output
end



# DON'T MODIFY THIS FILE -- the source is in file General Utilities.ipynb. Look there for further documentation and examples of running the code.



"""
safe_axes(axh; further_params...)

If you're going to make axh the current axes, this function
first makes axh's figure the current figure. Some Julias
break without that.

Any optional keyword-value arguments are passed on to axes()
"""

function safe_axes(axh; further_params...)
    figure(axh[:figure][:number])
    axes(axh; Dict(further_params)...)
end


"""
    ax = axisWidthChange(factor; lock="c", ax=nothing)
"""
function axisWidthChange(factor; lock="c", ax=nothing)
    if ax==nothing; ax=gca(); end
    x, y, w, h = ax[:get_position]()[:bounds]

    if lock=="l";
    elseif lock=="c" || lock=="m"; x = x + w*(1-factor)/2;
    elseif lock=="r"; x = x + w*(1-factor);
    else error("I don't know lock type ", lock)
    end

    w = w*factor;
    ax[:set_position]([x, y, w, h])

    return ax
end


"""
ax = axisHeightChange(factor; lock="c", ax=nothing)
"""
function axisHeightChange(factor; lock="c", ax=nothing)
    if ax==nothing; ax=gca(); end
    x, y, w, h = ax[:get_position]()[:bounds]

    if lock=="b";
    elseif lock=="c" || lock=="m"; y = y + h*(1-factor)/2;
    elseif lock=="t"; y = y + h*(1-factor);
    else error("I don't know lock type ", lock)
    end

    h = h*factor;
    ax[:set_position]([x, y, w, h])

    return ax
end


"""
   ax = axisMove(xd, yd; ax=nothing)
"""
function axisMove(xd, yd; ax=nothing)
    if ax==nothing; ax=gca(); end
    x, y, w, h = ax[:get_position]()[:bounds]

    x += xd
    y += yd

    ax[:set_position]([x, y, w, h])
    return ax
end


"""
[] = remove_xtick_labels(ax=NaN)

Given an axis object, or an array of axes objects, replaces each xtick label string with the empty string "".

If no axis is passed, uses gca() to work with the current axis.


"""
function remove_xtick_labels(ax=nothing)

    if ax==nothing
        ax = gca()
    end

    if typeof(ax) <: Array
        for i=1:length(ax)
            remove_xtick_labels(ax[i])
        end
        return
    end

    nlabels = length(ax[:xaxis][:get_ticklabels]())

    newlabels = Array{String,1}(nlabels)
    for i=1:length(newlabels);
        newlabels[i] = ""
    end

    ax[:xaxis][:set_ticklabels](newlabels)
    return
end



"""
[] = remove_ytick_labels(ax=NaN)

Given an axis object, or an array of axes objects, replaces each ytick label string with the empty string "".

If no axis is passed, uses gca() to work with the current axis.


"""
function remove_ytick_labels(ax=nothing)

    if ax==nothing
        ax = gca()
    end

    if typeof(ax) <: Array
        for i=1:length(ax)
            remove_ytick_labels(ax[i])
        end
        return
    end

    nlabels = length(ax[:yaxis][:get_ticklabels]())

    newlabels = Array{String,1}(nlabels)
    for i=1:length(newlabels);
        newlabels[i] = ""
    end

    ax[:yaxis][:set_ticklabels](newlabels)
    return
end




# DON'T MODIFY THIS FILE -- the source is in file General Utilities.ipynb. Look there for further documentation and examples of running the code.


using PyCall
# To specifically test with QT or Tk backends, uncomment one of the two following lines:
# (the pygui call must come *before* the very first using PyPlot call)
#
# pygui(:qt)
# pygui(:tk)
#
# Note that get_current_fig_position() and set_current_fig_position() work only with either QT or Tk at this point.
#
using PyPlot


# If the Python path does not already have the local directory in it
if PyVector(pyimport("sys")."path")[1] != ""
    # Then the following line is PyCall-ese for "add the current directory to the Python path"
    pushfirst!(PyVector(pyimport("sys")."path"), "")
end
# We use Python to enable callbacks from the figures:
kbMonitorModule = pyimport("kbMonitorModule")



__permanent_BP_store = []   # The user doesn't need to worry about this variable, it is here to ensure that
                            # kbMonitorModule.kb_monitor objects created inside the install_nearest_callback()
                            # function do not get deleted upon exit of that function

"""
BP = install_nearest_point_callback(fighandle, user_callback; user_data=nothing)

This function makes the figure indicated by fighandle interactive: any time the mouse is clicked
while pointing inside any of the axes of the figure, the function user_callback() will be called,
and will be passed parameters that indicate which of the points drawn in the axis is the one
closest to the clicked point (closest in Euclidean data units).

**WARNING** because this runs through PyCall, any errors in your user_callback function will sadly
not show up.  The function will simply fail to work. So be careful and debug with lots of print statements.


# PARAMETERS:

- fighandle       A matplotlib figure handle, e.g. the result of figure(2)

- user_callback   A function, which must take 4 or 5 parameters (see below). These will be passed to it as:


* PARAMETERS OF YOUR FUNCTION USER_CALLBACK:

        - xy          A 2-element tuple, indicating the (x,y) position of the drawn point closest to the clicked point

        - r           A scalar, indicating the Euclidean distance between the clicked location and xy

        - linehandle  A matplotlib handle to a Lines2D object (e.g., the result of plot([1,2], [3, 10])) or a PathCollection object (as returned by scatter()).

        - axhandle    A matplotlib handle to the axis (e.g., the result of gca()) in which the event occurred.

        - user_data   If `install_nearest_point_callback()` was called user_data set to something, then
                      your function will be called with *five* parameters, and the last one will be the contents of
                      the user_data

# OPTIONAL PARAMETERS FOR INSTALL_NEAREST_POINT_CALLBACK():

- user_data       Data to be stored internally and made available to the callback function. Default
                  is nothing, in which case the callback function is called with 4 params


# RETURNS:

- BP     A PyCall.PyObject kbMonitorModule_kb_monitor object. This object contains the underlying engine
linking the figure to the callback function. To disconnect that link, call "remove_BP(BP)". To disconnect
*all* existing BP-function links, call "remove_all_BPs()".


# EXAMPLE:

pygui(true)

```jldoctest
function mycallback(xy, r, h, ax)
    @printf("(%.3f,%.3f), r=%3f ", xy[1], xy[2], r);
    print(h)
    print(ax)
    print("\n")
end

BP = install_nearest_point_callback(figure(2), mycallback)
plot([2,2])
```

"""
function install_nearest_point_callback(fighandle, user_callback; user_data=nothing)

    function point_nearest_to_click(BP)
        bpe = BP[:buttonlist]()
        # Remove any leading clicks that weren't inside an axis:
        while length(bpe)>0 && ((bpe[1][1]==nothing) || (bpe[1][1]==Void))
            bpe = bpe[2:end]
        end
        if length(bpe)>0
            ax = BP[:buttonlist]()[1][1]   # the axis we're working with
            x  = BP[:buttonlist]()[1][2]   # the x of the clicked point
            y  = BP[:buttonlist]()[1][3]   # the y of the clicked point

            ch = ax[:get_children]()       # all children of the axis

            idx    = nothing    # this'll be the index of the data point closest to the clickpoint
            minJ   = nothing    # the smallest squared distance between data point and clickpoint found so far
            handle = nothing    # the matplotlib handle of the line object for which the closes data point is found
            dx     = nothing    # closest data point x position
            dy     = nothing    # closest data point y position

            # Look over all children of the axis:
            for i=1:length(ch)
                # But only consider line objects:
                if contains(pystring(ch[i]), "lines.Line2D")
                    D = ch[i][:get_data]()    # D will be a Tuple with xdata, ydata vectors
                elseif contains(pystring(ch[i]), "PathCollection")
                    D = ch[i][:get_offsets]()    # D will be a matrix with xdata, ydata columns
                    D = (D[:,1], D[:,2])         # Turn it into a Tuple like for Line2D objects
                end
                if contains(pystring(ch[i]), "lines.Line2D") || contains(pystring(ch[i]), "PathCollection")
                    J = (D[1] - x).^2 + (D[2] - y).^2
                    ix = indmin(J)
                    if idx == nothing || J[ix] < minJ   # if we did not yet have a minimum candidate or this one is better
                        idx = ix; minJ = J[ix]; handle = ch[i]   # store our candidate
                        dx = D[1][ix]; dy = D[2][ix]
                    end
                end
            end

            # @printf("install: Am about to call the user callback\n")
            if minJ != nothing
                if BP[:get_userdata]() == nothing
                    user_callback((dx,dy), sqrt(minJ), handle, ax)
                else
                    user_callback((dx,dy), sqrt(minJ), handle, ax, BP[:get_userdata]())
                end
            end
            # @printf("install: Just returned from the user callback\n")

            # After dealing with all the buttonclick callbacks and so on, bring focus back to the figure that was clicked:
            figure(ax[:figure][:number])
        end

        # We've dealt with the buttonclick, clear the buttonlist
        # @printf("Am about to clear the button list on button "); print(BP); print("\n")
        BP[:clear_buttonlist]()
    end

    BP = kbMonitorModule.kb_monitor(fighandle, callback = point_nearest_to_click, userData=user_data)
    global __permanent_BP_store = [__permanent_BP_store ; BP]

    return BP
end


"""
    userdata = install_callback_reporter(xy, r, axhandle, dothandle, userdata)

Useful as a debugging tool for `install_nearest_point_callback()`: can be used
with that function as a callback; when called, simply prints out its parameters.
Returns userdata, does not print it.

"""
function install_callback_reporter(xy, r, linehandle, axhandle, userdata=nothing)
    @printf("xy=(%g,%g), r=%g\n", xy[1], xy[2], r)
    print("Line Handle:\n"); print(linehandle); print("\n")
    print("Axis Handle:\n"); print(axhandle); print("\n")
    # print("User Data:\n"); print(userdata); print("\n")

    return userdata
end



"""
    remove_BP(BP::PyCall.PyObject)

Disconnects a kbMonitorModule.kb_monitor object from its figure
"""
function remove_BP(BP::PyCall.PyObject)
    if contains(pystring(BP), "kbMonitorModule.kb_monitor")
        BP[:__del__]()

        i = find(__permanent_BP_store .== BP)
        if length(i)>0;
            i = i[1];
            global __permanent_BP_store = __permanent_BP_store[[1:(i-1) ; (i+1):end]]
        end
    end
end


"""
    remove_all_BPs()

    Disconnects *all* kbMonitorModule.kb_monitor objects from their figures.
"""
function remove_all_BPs()
    for BP in __permanent_BP_store
        BP[:__del__]()
    end

    global __permanent_BP_store = []
end



# DON'T MODIFY THIS FILE -- the source is in file General Utilities.ipynb. Look there for further documentation and examples of running the code.




"""

    (x, y, w, h) = get_current_fig_position()

Works only when pygui(true) and when the back end is QT. Has been tested only with PyPlot.
"""
function get_current_fig_position()
    # if !contains(pystring(plt[:get_current_fig_manager]()), "FigureManagerQT")
    try
        if contains(pystring(plt[:get_current_fig_manager]()), "Tk")
            g = split(plt[:get_current_fig_manager]()[:window][:geometry](), ['x', '+'])
            w = parse(Int64, g[1])
            h = parse(Int64, g[2])
            x = parse(Int64, g[3])
            y = parse(Int64, g[4])
        elseif contains(pystring(plt[:get_current_fig_manager]()), "QT")
            x = plt[:get_current_fig_manager]()[:window][:pos]()[:x]()
            y = plt[:get_current_fig_manager]()[:window][:pos]()[:y]()
            w = plt[:get_current_fig_manager]()[:window][:width]()
            h = plt[:get_current_fig_manager]()[:window][:height]()
        else
            error("Only know how to work with matplotlib graphics backends that are either Tk or QT")
        end

        return (x, y, w, h)
    catch
        error("Failed to get current figure position. Is pygui(false) or are you using a back end other than QT or Tk?")
    end
end

"""

    set_current_fig_position(x, y, w, h)

Works only when pygui(true) and when the back end is QT. Has been tested only with PyPlot.
"""
function set_current_fig_position(x, y, w, h)
    # if !contains(pystring(plt[:get_current_fig_manager]()), "FigureManagerQT")
    try
        if contains(pystring(plt[:get_current_fig_manager]()), "Tk")
            plt[:get_current_fig_manager]()[:window][:geometry](@sprintf("%dx%d+%d+%d", w, h, x, y))
        elseif contains(pystring(plt[:get_current_fig_manager]()), "QT")
            plt[:get_current_fig_manager]()[:window][:setGeometry](x, y, w, h)
        else
            error("Only know how to work with matplotlib graphics backends that are either Tk or QT")
        end
    catch
        error("Failed to set current figure position. Is pygui(false) or are you using a back end other than QT?")
    end
end


"""
    C = capture_current_figure_configuration()

Collects the positions of all current figures and
prints out to the screen code, that can be copy-pasted,
that would reproduce that positioning configuration.

# PARAMETERS:

None

# RETURNS:

- C    A matrix that is nfigures-by-5 in size. You probably
    don't want this, you probably want the text printed to
    the screen, but here just in case.  Each row will have,
    in order: figure number, x, y, width, height
"""
function capture_current_figure_configuration()
    @printf("The following code will reproduce your current figure placement:\n\n")
    C = []
    for f in sort(plt[:get_fignums]())
        figure(f)
        x, y, w, h = get_current_fig_position()
        @printf("figure(%d); set_current_fig_position(%d, %d, %d, %d)   # x, y, width, height\n",
            f, x, y, w, h)
        C = [C ; [f x y w h]]
    end
    return C
end


end   # ====== END MODULE =============================
