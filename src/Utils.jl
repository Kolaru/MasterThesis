module Utils

using IntervalArithmetic
using IntervalRootFinding

import JSON

export UNITINTERVAL
export lower

const UNITINTERVAL = 0..1

function JSON.lower(rt::Root)
    X = interval(rt)
    N = length(X)
    bounds = zeros(2, N)
    for (k, x) in enumerate(X)
        bounds[1, k] = x.lo
        bounds[2, k] = x.hi
    end
    return Dict("bounds" => bounds, "status" => rt.status)
end

end
