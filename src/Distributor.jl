function PowerDistribute!(f::Field{Dom}, ps::Field{PowerDomain{Dom}}, out::Field{Dom}) where Dom <: StructuredDomain#{harmonic = true}
    pindices = get_pindices(ps)#.domain._pindices
    for i in eachindex(f)
        #out[i] = ps.val[pindices[i]]*f.val[i]
        f.val[i] *= ps.val[pindices[i]]
        #out[i] = f[i]
    end
    return nothing
end

function PowerDistribute!(f::Field{Dom}, ps::Field{PowerDomain{Dom}}) where Dom <: StructuredDomain
    PowerDistribute(f, ps, f)
    return nothing
end

function PowerDistribute(f::Field{Dom}, ps::Field{PowerDomain{Dom}}) where Dom
    pindices = ps.domain._pindices
    out = field(f.domain, similar(f.val))
    for i in eachindex(f)
        #out.val[i] = ps.val[pindices[i]]*f.val[i]
        out.val[i] = ps.val[pindices[i]]
    end
end
