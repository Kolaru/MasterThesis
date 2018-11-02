module Utils

using IntervalArithmetic
using IntervalRootFinding

import JSON

export UNITINTERVAL
export lower

const UNITINTERVAL = 0..1

function JSON.lower(rt::Root)
    return Dict("bounds" => JSON.lower(rt.interval), "status" => rt.status)
end

JSON.lower(X::Interval) = [X.lo, X.hi]

function JSON.lower(X::IntervalBox{N, T}) where {N, T}
    bounds = zeros(2, N)
    for (k, x) in enumerate(X)
        bounds[1, k] = x.lo
        bounds[2, k] = x.hi
    end
    return bounds
end

end
