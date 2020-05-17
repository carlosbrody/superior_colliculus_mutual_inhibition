using Optim
using TanhWalls

function f(x)
   return sum(x.^2)
end

bfunc, old2new, new2old = boundedVersion(f, [1 3; -100 100])

g = x -> ForwardDiff.gradient(bfunc, x)

result = optimize(
    bfunc, g, # h, seems overall faster without?
    old2new([2.0,2]), # NewtonTrustRegion();  seems overall faster without?
    Optim.Options(store_trace=true, show_trace=true, time_limit=500);
    inplace=false);
