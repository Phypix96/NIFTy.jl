#TODO macro for cleaner type signature
function PowerDistribute!(f::Field{T, N, S, Dom},
                          ps::Field{Tp, Np, Sp, PowerDomain{N, Dom}},
                          out::Field{T, N, S, Dom}) where {T, N, S, Dom <: StructuredDomain, Tp, Np, Sp}
    pind = pindices(ps)
    for i in eachindex(f)
        out[i] = f[i]*ps[pind[i]]
    end
    return nothing
end

function PowerDistribute!(f::Field{T, N, S, Dom}, ps::Field{Tp, Np, Sp, PowerDomain{N, Dom}}) where {T, N, S, Dom <: StructuredDomain, Tp, Np, Sp}
    PowerDistribute(f, ps, f)
    return nothing
end

function PowerDistribute(f::Field{T, N, S, Dom}, ps::Field{Tp, Np, Sp, PowerDomain{N, Dom}}) where {T, N, S, Dom, Tp, Np, Sp}
    out = similar(f)
    PowerDistribute!(f, ps, out)
    return out
end
