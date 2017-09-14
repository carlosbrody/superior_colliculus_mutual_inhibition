
"""
function x, cost, iters_used, last_Delta_x = one_d_minimizer(seed, func; tol=1e-5, maxiter=100, start_eta=0.1)

Minimizes a 1-d function using constrained Hessian minimization. 
We don't trust the long-range info from the Hessian too much, meaning that there's
a given (adaptive) step size. If Newton's method suggests a step smaller than that step
size, we take it. Otherwise, we only move in the direction of the gradient by the stepsize.

Adaptive step size: Every step that the cost function gets smaller, the step grows by a factor 
of 1.1. If a step would have led to a larger cost function, the step is not taken, and
the size falls by a factor of 2.

PARAMETERS:
===========

seed       A float, the starting point for the minimization

func       A one dimensional function, takes a float and returns a float that represents the current cost.


OPTIONAL PARAMETERS:
====================

tol=1e-5   If a step would lead to a change in cost that is smaller in magnitude than tol, stop the minimization.

maxiter=100   Maximum number of iterations for the minimization

start_eta=0.1  The starting value of the step size


RETURNS:
========

x0     the value of x that minimizes f(x)

cost   The value of f(x0) 

niters  The number of iterations done

Dparam  The last change to the parameter value

"""
function one_d_minimizer(seed, func; tol=1e-5, maxiter=100, start_eta=0.1)
    eta = start_eta;
    lambdavec = [seed]


    out = DiffBase.HessianResult(lambdavec)
    ForwardDiff.hessian!(out, func, lambdavec)
    cost = DiffBase.value(out)
    grad = DiffBase.gradient(out)
    hess = DiffBase.hessian(out)

    i = delta_lambda = 0;  # declare it so we can use it after the if
    for i in [1:maxiter;]
        h_delta = - grad./hess;
        if abs(h_delta[1]) < eta
            new_lambdavec = lambdavec + h_delta
        else
            new_lambdavec = lambdavec - eta*sign(grad)
        end

        delta_lambda = new_lambdavec[1] - lambdavec[1]
        if abs(delta_lambda) < tol
            break
        end
        
        ForwardDiff.hessian!(out, func, new_lambdavec)
        new_cost = DiffBase.value(out)
        new_grad = DiffBase.gradient(out)
        new_hess = DiffBase.hessian(out)

        if new_cost .< cost
            lambdavec[1] = new_lambdavec[1]
            cost = new_cost;
            grad = new_grad;
            hess = new_hess;

            eta = eta*1.1
        else
            eta = eta/2
        end

        # @printf "%d: cost=%.3f lambda=%.3f\n" i cost lambdavec[1]
    end
    
    return lambdavec[1], cost, i, delta_lambda
end
    
