####################################################################################################
####################################################################################################
#Simple Pointwise Operators
exp(dom::AbstractDomain) = PointwiseOperator(dom, [exp], [exp])
expm1(dom::AbstractDomain) = PointwiseOperator(dom, [expm1], [exp])
log(dom) = PointwiseOperator(dom, [log], [inv])
log1p(dom) = PointwiseOperator(dom, [log1p], [x -> 1/(1 + x)])
sin(dom) = PointwiseOperator(dom, [cos], [sin])
cos(dom) = PointwiseOperator(dom, [cos], [sin, -])

add(dom, val::Union{Number, AbstractArray{Number}}) = PointwiseOperator(dom, [x -> +(x, val)], [identity])
scale(dom, val::Union{Number, AbstractArray{Number}}) = PointwiseOperator(dom, [x -> *(x, val)], [x -> *(x, val)])

####################################################################################################
####################################################################################################
#Simple Linear Operators
struct HartleyOperator <: LinearOperator
    domain::StructuredDomain{N, true} where N
end
target(op::HartleyOperator) = getcodomain(op.domain)
apply(op::HartleyOperator, val) = hartley(val)
adjoint(op::HartleyOperator, val) = apply(op, val)

hartley(op::AbstractOperator) = combine_operators(HartleyOperator(target(op)), op)


struct FourierOperator <: LinearOperator
    domain::StructuredDomain{N, true} where N
end
target(op::FourierOperator) = getcodomain(op.domain)
apply(op::FourierOperator, val) = fft(val)
#adjoint(op::FourierOperator, val) = fft(conj.(val))

####################################################################################################
####################################################################################################
#Multi-Field Operators
struct AdditionOperator <: AbstractOperator
    domain
    target
    operators
end

struct SubtractionOperator <: AbstractOperator
    domain
    target
    operators
end

struct MultiplicationOperator <: AbstractOperator
    domain
    target
    operators
end

struct DivisionOperator <: AbstractOperator
    domain
    target
    operators
end

function +(op1::AbstractOperator, op2::AbstractOperator)
    @assert target(op1) == target(op2)
    dom = (domain(op1), domain(op2))
    return AdditionOperator(dom, target(op1), (op1, op2))
end

function -(op1::AbstractOperator, op2::AbstractOperator)
    @assert target(op1) == target(op2)
    dom = (domain(op1), domain(op2))
    return SubtractionOperator(dom, target(op1), (op1, op2))
end

function *(op1::AbstractOperator, op2::AbstractOperator)
    @assert target(op1) == target(op2)
    dom = (domain(op1), domain(op2))
    return MultiplicationOperator(dom, target(op1), (op1, op2))
end

function /(op1::AbstractOperator, op2::AbstractOperator)
    @assert target(op1) == target(op2)
    dom = (domain(op1), domain(op2))
    return DivisionOperator(dom, target(op1), (op1, op2))
end


function apply(op::AdditionOperator, vals)
    #TODO replace zeros(size(target)) with zeros(target)
    res = zeros(size(target(op)))
    for (operation, val) in zip(op.operators, vals)
        res .+= apply(operation, val)
    end
    return res
end

function apply(op::SubtractionOperator, (val1, val2))
    return apply(op.operators[1], val1) .- apply(op.operators[2], val2)
end

function apply(op::MultiplicationOperator, vals)
    #TODO replace zeros(size(target)) with zeros(target)
    res = ones(size(target(op)))
    for (operation, val) in zip(op.operators, vals)
        res .*= apply(operation, val)
    end
    return res
end

function apply(op::SubtractionOperator, (val1, val2))
    return apply(op.operators[1], val1) ./ apply(op.operators[2], val2)
end
end
####################################################################################################
####################################################################################################
#Convenience Methods
for f in [:exp, :expm1, :log, :log1p, :sin, :cos]
    @eval $f(op::AbstractOperator) = combine_operators($f(target(op)), op)
end

