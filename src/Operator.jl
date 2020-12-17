#TODO handle identity separately in pointwise operators


abstract type AbstractOperator end

struct Operator <: AbstractOperator
    domain::AbstractDomain
    target::AbstractDomain
    operation
    gradient
    reverse
end

struct LinearOperator <: AbstractOperator
    domain::AbstractDomain
    target::AbstractDomain
    operations::Array
    adjoints::Array
end

struct PointwiseOperator <: AbstractOperator
    domain::AbstractDomain
    operations::Array
    gradients::Array
end

struct OperatorChain <: AbstractOperator
    domain::Union{AbstractDomain, Tuple{Vararg{AbstractDomain}}}
    target::Union{AbstractDomain, Tuple{Vararg{AbstractDomain}}}
    operators::Array{O} where O <: AbstractOperator
end


####################################################################################################
####################################################################################################
#Property interfaces
domain(op::Operator) = op.domain
target(op::Operator) = op.target

domain(op::LinearOperator) = op.domain
target(op::LinearOperator) = op.target

domain(op::PointwiseOperator) = op.domain
target(op::PointwiseOperator) = op.domain

domain(op::OperatorChain) = op.domain
target(op::OperatorChain) = op.target

get_operators(op::AbstractOperator) = identity(op)
get_operators(op::OperatorChain) = op.operators

####################################################################################################
####################################################################################################
#combine_operators Interface
function (op_later::AbstractOperator)(op_first::AbstractOperator)
    return combine_operators(op_later, op_first)
end


function combine_operators(op_later::LinearOperator, op_first::LinearOperator)
    @assert target(op_first) == domain(op_later)
    operations = vcat(op_first.operations, op_later.operations)
    adjoints = vcat(op_later.adjoints, op_first.adjoints)
    return LinearOperator(op_later.domain, op_first.target, operations, adjoints)
end

function combine_operators(op_later::PointwiseOperator, op_first::PointwiseOperator)
    @assert domain(op_later) == domain(op_first)
    operations = vcat(op_first.operations, op_later.operations)
    gradients = vcat(op_first.gradients, op_later.gradients)
    return PointwiseOperator(op_first.domain, operations, gradients)
end

function combine_operators(op_later, op_first)
    @assert target(op_first) == domain(op_later)
    operators = vcat(get_operators(op_first), get_operators(op_later))
    return OperatorChain(domain(op_first), target(op_later), operators)
end


####################################################################################################
####################################################################################################
#apply Interface

#TODO generalize to NTuple for val
function (op::AbstractOperator)(val)
    res = deepcopy(val)
    return apply!(op, val)
end

function (op::AbstractOperator)(vals...)
    res = deepcopy(vals)
    res = apply!(op, res)
    return res
end

function apply!(op::LinearOperator, val::AbstractArray)
    for operation in op.operations
        val = operation(val)
    end
    return val
end

#TODO generalize to NTuple for val
function apply!(op::PointwiseOperator, val::AbstractArray)
    for operation in op.operations
        @avx @. val = operation(val)
    end
    return val
end

function apply!(op::OperatorChain, val)
    for operator in op.operators
        val = apply!(operator, val)
    end
    return val
end


#function apply!(op::AdditionOperator, val1, val2)
#    @avx @. val1 += val2
#end
#function apply!(op::MultiplicationOperator, val1, val2)
#    @avx @. val1 *= val2
#end
####################################################################################################
####################################################################################################
#apply linearization Interface

#function (op::Operator)(val, jac)
#    new_jac(x) = op.jacobian(val, jac(x))
#    new_val = op.operation(val)
#    return (new_val, new_jac)
#end
#
#function (op::LinearOperator)(val, jac)
#    res = op(val)
#    jacobian(x) = op(jac(x))
#    return res, jacobian
#end
#
#function (op::PointwiseOperator)(val, jac)
#    gradient = ones(size(val))
#    res = copy(val)
#    for (operation, grad) in zip(op.operations, op.gradients)
#        @avx @. gradient *= grad(res)
#        @avx @. res = operation(res)
#    end
#    jacobian(x) = @avx gradient .* jac(x)
#    return res, jacobian
#end
#
#function (op::OperatorChain)(val, jac)
#    res = copy(val)
#    for operator in op.operators
#        #TODO here, res can be modified by apply, so pointwise operators don't need to make a copy
#        res, jac = operator(res, jac)
#    end
#    return res, jac
#end

####################################################################################################
####################################################################################################
#miscellaneous

adjoint(op::LinearOperator) = LinearOperator(target(op), domain(op), op.adjoints, op.operations)



exp(dom::AbstractDomain) = PointwiseOperator(dom, [exp], [exp])
exp(op::AbstractOperator) = PointwiseOperator(target(op), [exp], [exp])(op)
log(dom) = PointwiseOperator(dom, [log], [inv])
sin(dom) = PointwiseOperator(dom, [cos], [sin])
cos(dom) = PointwiseOperator(dom, [cos], [-, sin])
add(dom, val::Number) = PointwiseOperator(dom, [x -> +(x, val)], [identity])
scale(dom, val::Number) = PointwiseOperator(dom, [x -> *(x, val)], [identity])


hartley(dom) = LinearOperator(dom, getcodomain(dom), [hartley], [hartley])
