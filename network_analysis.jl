# Build the connectivity matrix from the parameter vector. Can optionally include an extra leak on the diagonal
function build_w(params; ofs=zeros(4,1))
    w = [params[6]-ofs[1]   params[16]          params[10]          params[11];
         params[3]          params[9]-ofs[2]    params[5]           params[13];
         params[13]         params[5]           params[9]-ofs[3]    params[3];
         params[11]         params[10]          params[16]          params[6]-ofs[4]];
    return w;
end

# gets the schur decomposition for a model solution by loading the file, and making the connectivity matrix
function get_schur(filename;ofs=zeros(4,1))
    # get parameters
    f = load(filename);    

    # build W
    w = build_w(f["pars3"];ofs=ofs);

    # get schur
    s = schur(w);
    return s;
end

# test a random matrix for its schur properties. Can either be iid randn, or impose symmetry, or impose symmetry with parameter bounds
function test_random(num2test;add_symmetry=false, add_bounds=false)
numcomplex = 0;
full_set = 0;
counts = zeros(5,2);
for i=1:num2test
    if add_symmetry & add_bounds
        params = rand(16,1).*6 - 3;
        params[6] = abs(params[6]);
        params[9] = abs(params[9]);
        # all are -3,3 ,except self weights, which are 0 to 3
        w = build_w(params);
    elseif add_symmetry
        params = randn(16,1);
        w = build_w(params);
    else
        w = randn(4,4)
    end
    s = schur(w);
    if !isreal(s[3]);   numcomplex +=1; end
    dexes = [];
    for j=1:4        
        mode = s[2][:,j];
        dex =0;
        if (all( (mode .> 0) .== [false;false;false;false])) | (all( (mode .> 0) .== [true;true;true;true]))
            dex = 1; # All active mode
        elseif (all( (mode .> 0) .== [true;true;false;false])) | (all( (mode .> 0) .== [false;false;true;true]))
            dex = 2; # side mode
        elseif (all( (mode .> 0) .== [true;false;false;true])) | (all( (mode .> 0) .== [false;true;true;false]))
            dex = 3; # pro/anti mode
        elseif (all( (mode .> 0) .== [true;false;true;false])) | (all( (mode .> 0) .== [false;true;false;true])) 
            dex = 4; # diagonal mode
        else 
            dex = 5;
        end
        if real(s[3][j]) > 0
            stable = 1;
        else
            stable = 2;
        end
        counts[dex, stable] +=1;
        dexes = [dexes... dex];
        if length(unique(dexes)) == 4
            full_set +=1;
        end    
    end
end
return full_set, numcomplex, counts
end


# get schur decomposition for every solution
function tabulate_modes(results;ofs=zeros(4,1))
    schurs = [];
    complex_list = [];
    counts = zeros(5,3);
    full_set = [];
    interactions = zeros(4,4,2);

    for i=1:length(results["files"])
        s = get_schur(results["files"][i];ofs=ofs);
        schurs= [schurs..., s];
        if !isreal(s[3])
           complex_list = [complex_list..., true];
        else
            complex_list = [complex_list..., false];
        end
        dexes = [];
        for j=1:4        
            mode = s[2][:,j];
            dex =0;
            if (all( (mode .> 0) .== [false;false;false;false])) | (all( (mode .> 0) .== [true;true;true;true]))
                dex = 1; # All active mode
            elseif (all( (mode .> 0) .== [true;true;false;false])) | (all( (mode .> 0) .== [false;false;true;true]))
                dex = 2; # side mode
            elseif (all( (mode .> 0) .== [true;false;false;true])) | (all( (mode .> 0) .== [false;true;true;false]))
                dex = 3; # pro/anti mode
            elseif (all( (mode .> 0) .== [true;false;true;false])) | (all( (mode .> 0) .== [false;true;false;true])) 
                dex = 4; # diagonal mode
            else 
                dex = 5;
            end
     #       if !isreal(s[3][j])
     #           stable = 3;
            if real(s[3][j]) > 0
                stable = 1;
            else
                stable = 2;
            end
            counts[dex, stable] +=1;
            dexes = [dexes... dex];
        end

        if length(unique(dexes)) == 4
            full_set =[full_set..., true];
        else
            full_set =[full_set..., false];       
        end

if false    
        # Analyze interactions    
        if length(unique(dexes)) == 4
            for j=1:4
                for k=1:4
                    pos = s[1][j,k] > 1e-4;
                    neg = s[1][j,k] < -1e-4;
                    if pos
                        interactions[dexes[j],dexes[k],1] += 1;
                    elseif neg
                        interactions[dexes[j],dexes[k],2] += 1;
                    end
                end
            end
        end
end
    end

    return full_set, complex_list, counts, schurs, interactions
end

# Calculate statistics for different null models
random_full_set, random_numcomplex, random_counts = test_random(10000;add_symmetry=false, add_bounds=false);
sym_full_set,    sym_numcomplex,    sym_counts    = test_random(10000;add_symmetry=true,  add_bounds=false);
bound_full_set,  bound_numcomplex,  bound_counts  = test_random(10000;add_symmetry=true,  add_bounds=true);

# get the schur results for all the model solutions
results = load_farm_cost_filter("C32", "MiniC32"; threshold = -0.0001)
full_set, numcomplex, counts, schurs, interactions = tabulate_modes(results;ofs=[0;0;0;0]);

