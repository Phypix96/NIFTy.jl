function PowerDistribute!(f::Field{Dom}, ps::Field{PowerDomain{Dom}}, out::Field{Dom}) where Dom <: StructuredDomain#{harmonic = true}
    pindices = ps.domain._pindices
    for i in eachindex(f)
        out[i] = ps.val[pindices[i]]*f.val[i]
    end
    return nothing
end

function PowerDistribute!(f::Field{Dom}, ps::Field{PowerDomain{Dom}}) where Dom <: StructuredDomain
    PowerDistribute(f, ps, f)
    return nothing
end
