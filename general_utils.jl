# DON'T MODIFY THIS FILE -- the source is in file General Utilities.ipynb. Look there for further documentation and examples of running the code.



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



"""
function print_vector(vec)

Takes a vector and uses @printf to put it on the screen with [%.3f, %.3f] format. 

If passed a symbol (which must evaluate to a vector), then prints the string for that symbol,
an equals sign, the vector, and ends by adding a carriage return \n.
"""
function print_vector(vec)

    if typeof(vec)==Symbol
        mystr = string(vec)
        @printf("%s = ", mystr);
        print_vector(eval(vec))
        @printf("\n");
        return
    end
    
    @printf "["
    for p in [1:length(vec);]
        @printf("%.3f", vec[p])
        if p < length(vec) @printf ", "; end
    end
    @printf "]"
end


"""
function print_vector_g(vec)

Takes a vector and uses @printf to put it on the screen with [%g, %g] format. 

If passed a symbol (which must evaluate to a vector), then prints the string for that symbol,
an equals sign, the vector, and ends by adding a carriage return \n.
"""
function print_vector_g(vec)

    if typeof(vec)==Symbol
        mystr = string(vec)
        @printf("%s = ", mystr);
        print_vector_g(eval(vec))
        @printf("\n");
        return
    end
    
    @printf "["
    for p in [1:length(vec);]
        @printf("%g", vec[p])
        if p < length(vec) @printf ", "; end
    end
    @printf "]"
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



