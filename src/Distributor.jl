struct DistributionOperator <: AbstractOperator
    domain
    target
    pindices
end

domain(op::DistributionOperator) = op.domain
target(op::DistributionOperator) = op.target

function apply!(op::DistributionOperator, (val_full, val_distribute))
    for i in eachindex(val_full)
        val_full[i] *= val_distribute[op.pindices[i]]
    end
    return val_full
end


PowerDistributor(pdom::PowerDomain) = begin
    DistributionOperator((codomain(pdom), pdom), codomain(pdom), pdom._pindices)
end
