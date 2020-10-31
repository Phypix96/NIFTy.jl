#TODO macro for cleaner type signature
function PowerDistribute!(f::Field{T, N, Dom},
                          ps::Field{Tp, Np, PowerDomain{N, Dom}},
                          out::Field{T, N, Dom}) where {T, N, Dom <: StructuredDomain, Tp, Np}
    pind = pindices(ps)
    for i in eachindex(f)
        out[i] = f[i]*ps[pind[i]]
    end
    return nothing
end

function PowerDistribute!(f::Field{T, N, Dom}, ps::Field{Tp, Np, PowerDomain{N, Dom}}) where {T, N, Dom <: StructuredDomain, Tp, Np}
    PowerDistribute(f, ps, f)
    return nothing
end

function PowerDistribute(f::Field{T, N, Dom}, ps::Field{Tp, Np, PowerDomain{N, Dom}}) where {T, N, Dom, Tp, Np}
    out = similar(f)
    PowerDistribute!(f, ps, out)
    return out
end
