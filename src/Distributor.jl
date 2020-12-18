struct DistributionOperator <: AbstractOperator
    domain
    target
    pindices
end

domain(op::DistributionOperator) = op.domain
target(op::DistributionOperator) = op.target

function apply(op::DistributionOperator, (val_full, val_distribute))
    res = similar(val_full)
    @avx for i in eachindex(val_full)
        res[i] = val_full[i] * val_distribute[op.pindices[i]]
    end
    return res
end


PowerDistributor(pdom::PowerDomain) = begin
    DistributionOperator((codomain(pdom), pdom), codomain(pdom), pdom._pindices)
end
