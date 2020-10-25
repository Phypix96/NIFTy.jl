function PowerDistribute!(f::Field{sD}, ps::Field{PowerDomain}, out::Field{sD}) where sD <: StructuredDomain
    pindices = ps.domain._pindices
    for i in eachindex(f)
        out[i] = ps.val[pindices[i]]
    end
    return nothing
end

function PowerDistribute!(f::Field{sD}, ps::Field{PowerDomain})
    PowerDistribute(f, ps, f)
    return nothing
end
