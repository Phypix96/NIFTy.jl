#TODO macro for cleaner type signature
function PowerDistribute(f::Field{T, N, Dom},
                         ps::Field{Tp, Np, PowerDomain{N, Dom}}) where {T, N, Dom <: StructuredDomain, Tp, Np}
    pind = pindices(ps)
    out = similar(f)
    for i in eachindex(f)
        out[i] = f[i]*ps[pind[i]]
    end
    return out
end

