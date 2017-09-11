"""
[] = remove_xtick_labels(ax)

Given an axis object, or an array of axes objects, replaces each xtick label string with ""


"""
function remove_xtick_labels(ax)

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
new_fname = next_file(basename, ndigits)

Returns a numbered and presumably unused filename starting with the string basename, followed by an integer
digit. The returned integer will be one higher than the number of existing filenames starting with basename,
and will be written with ndigits numbers, using leading eros if necessary.

# EXAMPLE:

If there are already 8 files starting with "model" then

> next_file("model_", 4)

"model_0009"
"""
function next_file(basename, ndigits)
    fnames = readdir()
    matched_filenames = Array{Bool}(length(fnames))
    for i=1:length(fnames)
        matched_filenames[i] = ismatch(Regex(@sprintf("^%s", basename)), fnames[i])
    end
    
    mynum = length(find(matched_filenames))+1
    myname = @sprintf("%d", mynum)
    while length(myname)<ndigits
        myname = "0" * myname
    end

    return basename * myname
end
