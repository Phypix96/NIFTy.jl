#TODO handle identity separately in pointwise operators
abstract type AbstractOperator end

abstract type  LinearOperator <: AbstractOperator end

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
domain(op::AbstractOperator) = op.domain
target(op::AbstractOperator) = op.target

target(op::PointwiseOperator) = op.domain

get_operators(op::AbstractOperator) = identity(op)
get_operators(op::OperatorChain) = op.operators

####################################################################################################
####################################################################################################
#combine_operators Interface
function (op_later::AbstractOperator)(op_first::AbstractOperator)
    return combine_operators(op_later, op_first)
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

function (op::AbstractOperator)(f::Union{Field, MultiField})
    @assert f.domain == domain(op)
    return apply(op, f.val)
end

#function (op::AbstractOperator)(vals...)
#    return apply(op, vals)
#end

function apply(op::AbstractOperator, vals::NamedTuple)
    vals = [vals[key] for key in keys(domain(op))]
    return apply(op, vals)
end

#TODO generalize to NTuple for val
function apply(op::PointwiseOperator, val::AbstractArray)
    res = @avx @. op.operations[1](val)
    for operation in op.operations[2:end]
        @avx @. res = operation(res)
    end
    return res
end

function apply(op::OperatorChain, val)
    for operator in op.operators
        val = apply(operator, val)
    end
    return val
end


####################################################################################################
####################################################################################################
#apply linearization Interface

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
