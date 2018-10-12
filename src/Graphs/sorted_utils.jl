#=
This file contains utility function optimized for sorted sequences.
=#

"""
    in_sorted(vec::Vector{T}, x::T) where T

Return if `x` is in the sorted sequence `vec`.
"""
function in_sorted(vec::Vector{T}, x::T) where T
    searchsortedfirst(vec, x) > length(vec)
end

"""
    insert_sorted!(vec::Vector{T}, x::T) where T

Insert `x` in the sorted sequence `vec` such that it remains sorted.
"""
function insert_sorted!(vec::Vector{T}, x::T) where T
    insert!(vec, searchsortedfirst(vec, x), x)
end

"""
    replace_sorted!(vec::Vector{T}, xold::T, xnew::T) where T

Replace all occurrence of `xold` by `xnew` in the sorted sequence `vec`.
"""
function replace_sorted!(vec::Vector{T}, xold::T, xnew::T) where T
    sp = splice!(vec, searchsorted(vec, xold))
    sf = searchsortedfirst(vec, xnew)
    for i in 1:length(sp)
        insert!(v, sf, xnew)
    end
end

"""
    remove_sorted!(vec::Vector{T}, x::T) where T

Remove one occurence of `x` in the sorted sequecnce `vec`.
"""
function remove_sorted!(vec::Vector{T}, x::T) where T
    splice!(vec, searchsortedfirst(vec, x))
end

"""
    removeall_sorted!(vec::Vector{T}, x::T) where T

Remove all occurences of `x` in the sorted sequecnce `vec`.
"""
function removeall_sorted!(vec::Vector{T}, x::T) where T
    splice!(vec, searchsorted(vec, x))
end
