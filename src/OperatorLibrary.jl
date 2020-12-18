####################################################################################################
####################################################################################################
#Simple Pointwise Operators
exp(dom::AbstractDomain) = PointwiseOperator(dom, [exp], [exp])
expm1(dom::AbstractDomain) = PointwiseOperator(dom, [expm1], [exp])
log(dom) = PointwiseOperator(dom, [log], [inv])
log1p(dom) = PointwiseOperator(dom, [log1p], [x -> 1/(1 + x)])
sin(dom) = PointwiseOperator(dom, [cos], [sin])
cos(dom) = PointwiseOperator(dom, [cos], [sin, -])

add(dom, val::Number) = PointwiseOperator(dom, [x -> +(x, val)], [identity])
scale(dom, val::Number) = PointwiseOperator(dom, [x -> *(x, val)], [identity])

####################################################################################################
####################################################################################################
#Simple Linear Operators
hartley(dom) = LinearOperator(dom, getcodomain(dom), [hartley], [hartley])
fourier(dom) = LinearOperator(dom, getcodomain(dom), [fft], [fft])


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
        res .+= operation(val)
    end
    return res
end

function apply(op::SubtractionOperator, (val1, val2))
    return op.operators[1](val1) .- op.operators[2](val2)
end

function apply(op::MultiplicationOperator, vals)
    #TODO replace zeros(size(target)) with zeros(target)
    res = ones(size(target(op)))
    for (operation, val) in zip(op.operators, vals)
        res .*= operation(val)
    end
    return res
end

function apply(op::SubtractionOperator, (val1, val2))
    return op.operators[1](val1) ./ op.operators[2](val2)
end
####################################################################################################
####################################################################################################
#Convenience Methods
for f in [:exp, :expm1, :log, :log1p, :sin, :cos, :hartley, :fourier]
    @eval $f(op::AbstractOperator) = combine_operators($f(target(op)), op)
end

