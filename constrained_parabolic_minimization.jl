
"""
function constrained_parabolic_minimization(H, G, r; tol=1e-6, min_only=true, do_plot=false, 
    verbose=false, efactor=3.0, max_efactor_tries=10, 
    lambdastepsize=0.003, minimum_tol=1e-24, tol_delta=1e-3, maxiter=200)

Given a Hessian matrix, a gradient vector, and a desired radius from the origin, finds the vector 
that minimizes the parabola defined by the Hessian and the gradient, subject to the constraint that the
vector's length equals the desired radius.

# PARAMETERS:

- H      A square symmetric matrix. 

- G      A vector, length equal to the size(H,2)

- r      desired radius

# OPTIONAL PARAMETERS:

- tol=1e-6        Numerical tolerance on the computations

- min_only=true   Return only the minimum, or, if false, all xs for all lambdas that match x'*x = r^2 
                    (some of thos might be maxima, not minima)

- efactor=3       The initial exploration of lambdas will go from -efactor(max(absolute(eig(H)))) to +efactor(max(absolute(eig(H))))
                If that does not produce a solution within the indicates tolerance, then efactor is multiplied by efactor_growth
                and we try again, up to max_efactor_tries.

- efactor_growth  see efactor

- max_efactor_tries   see efactor

- minimum_tol     Each candidate value of lambda, identified as a minimum in the grid scan and 
                proposed as satisfying |x|^2=r^2, is then used to seed a Newton's method one-d minimization with a certain tolerance.
                If the results of that search do not satisfy |x|^2=r^2 within the desired tolerance, then the tolerance
                for the 1-d search is reduced (to make the one-d search more exact), by a factor of tol_delta, 
                up to minimu_tol

- tol_delta       see minimum_tol

- max_iter        Maximum number of iterations in the 1-d search (see minimum_tol).


- lambdastepsize=0.003    The step size for initial exploration of lambdas, un units of efactor. It sshould
                probably scale with the smallest difference in the eigenvalues of H; that has not been implemented yet.

- do_plot         If true, produces plot of the initial gridscans of lambda versus (x'*x - r^2)^2

- verbose         If true, produces diagnostic info as it goes.


RETURNS:
========

- x        The vector that minimizes 0.5*x'*H*x + x'*G subject to x'*x = r
 
- J        0.5*x'*H*x + x'*G at the returned x

- lambda   value of the Lagrange multiplier at which the radius constraint is satisfied

- c        The squared difference between the length of x and r. Should be small, otherwise somthing went wrong!

- niters   The number of iterations used in the one-d minimiation that identified the output lambdas.

- Dlambda  The last change in lambda during the one-d minimization. If the one-d minimization did not reach
         its maxiters, then this will be smaller than the one-d minimization's tolerance.

"""
function constrained_parabolic_minimization(H, G, r; tol=1e-6, min_only=true, 
    do_plot=false, verbose=false, efactor=3.0, efactor_growth=4, max_efactor_tries=10, 
    lambdastepsize=0.003, minimum_tol=1e-24, tol_delta=1e-3, maxiter=200)

    #  --- First a couple of helper functions ----
    
    """
    function x_of_lambda(lambda)

    Given square matrix H, vector G, and passed scalar lambda, returns the vector x that minimizes
    
    0.5 x'*H*x + x'*G - lambda *x'*x

    """
    function x_of_lambda(lambda)
        return inv(H - lambda*eye(size(H,1)))*(-G)
    end
    
    
    """
    function q(lambda, r)

    Returns the squared difference between r and the norm of x_of_lambda(lambda).
    """
    function q(lambda, r)
        return (r - norm(x_of_lambda(lambda)))^2
    end


    # efactor is the factor that multiplies the biggest eigenvalue of H, to determine the range over which we'll
    # look for a lambda that satisfies the norm(x)==r requirement. If we don't find a solution, we iteratively 
    # increase efactor to try to get there, for a maximum of max_efactor_tries
    for m=1:max_efactor_tries
        # First scan lambda to find good candidates for minimizing the parabolic 
        # surface under the x'*x = r^2 constraint
        L = eig(H)[1]
        L0 = maximum(abs(L))
        lambdas = L0*efactor*[-1.0:lambdastepsize:1.0;]
        costs = zeros(size(lambdas))
        for i in [1:length(lambdas);]
            try 
                costs[i] = q(lambdas[i], r)
            catch
                costs[i] = Inf
            end
        end

        if do_plot
            figure(2); clf();
            plot(lambdas, costs, "b.-")
            xlabel("lambda")
            ylabel("cost")
        end

        # Take all candidates where the derivative of costs changes sign 
        # from negative to positive (those would be minima),
        # plus the smallest and the largest lambdas tested, as candidates
        g = append!(prepend!(find(diff(sign(diff(costs))) .> 0.99), [1]), [length(lambdas)])
        if verbose
            @printf("cpm: g (candidate indices) are : ");           print_vector_g(g);        print("\n")
            @printf("cpm: and their corresponding costs are : ");   print_vector(costs[g]);   print("\n");
            @printf("cpm: and their corresponding lambdas are : "); print_vector(lambdas[g]); print("\n");
            @printf("cpm: the minimum cost was : %g\n", minimum(costs[g]))
        end
        # found_it_flag = 0  # A flag for when we've found at least one lambda that satisfies the r constraint
        mytol = tol

        while mytol > minimum_tol
            lambdas_out = zeros(size(g))
            costs_out   = zeros(size(g))
            niters_out  = zeros(size(g))
            Dlambda_out = zeros(size(g))
            for i in [1:length(g);]
                lambdas_out[i], costs_out[i], niters_out[i], Dlambda_out[i] = 
                one_d_minimizer(lambdas[g[i]], x -> q(x[1], r), start_eta=1, tol=mytol, maxiter=maxiter)
            end

            # Eliminate any lambdas where x'*x doesn't match our desired value r
            I = find(costs_out .< tol)
            lambdas_out = lambdas_out[I]; costs_out = costs_out[I];
            niters_out  = niters_out[I];  Dlambda_out = Dlambda_out[I]

            if length(I) > 0; break; end

            mytol *= tol_delta
        end
        if verbose
            @printf("%d : After searching for lambdas with efactor=%g, we found these : ", m, efactor)
            print_vector_g(lambdas_out); print("\n")
        end
        if length(lambdas_out) > 0; break; end;
        efactor = efactor*efactor_growth
    end
    
    # Eliminate any repeated lambdas, to within the specified numerical tolerance.
    I = setdiff(1:length(lambdas_out), find(diff(lambdas_out) .< tol))
    lambdas_out = lambdas_out[I]; costs_out = costs_out[I];
    niters_out  = niters_out[I];  Dlambda_out = Dlambda_out[I]
    
    # Find the parabolic estimate of the cost function at these points
    J  = zeros(size(lambdas_out))
    xs = zeros(length(G), length(lambdas_out))
    for i in [1:length(J);]
        xs[:,i] = x_of_lambda(lambdas_out[i])
        J[i] = (0.5*xs[:,i]'*H*xs[:,i] + xs[:,i]'*G)[1]
    end

    # Find and return only the x that has the smallest J
    if min_only
        I = indmin(J)    
    else
        I = 1:length(J)
    end
    return xs[:,I], J[I], lambdas_out[I], costs_out[I], niters_out[I], Dlambda_out[I]
end

