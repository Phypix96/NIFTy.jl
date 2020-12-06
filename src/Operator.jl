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

struct CombinedOperator <: AbstractOperator
    domain::AbstractDomain
    target::AbstractDomain
    operators::Array{O} where O <: AbstractOperator
end


function combine_operators(op1::LinearOperator, op2::LinearOperator)
    @assert op1.target == op2.domain
    operation(x) = op1.operation(op2.operation(x))
    adjoint(x) = op1.adjoint(op2.adjoint(x))
    return LinearOperator(op1.domain, op2.target, operation, adjoint)
end

function combine_operators(op1::PointwiseOperator, op2::PointwiseOperator)
    @assert op1.domain == op2.domain
    operation(x) = op1.operation(op2.operation(x))
    adjoint(x) = op2.adjoint(op1.adjoint(x))
    linearization(x) = begin 
        x, lin2 = op2.linearization(x)
        x, lin1= op1.linearization(x)
        jac(x) = lin1*lin2*x
        return x, jac
    end
    return PointwiseOperator(op1.domain, operation, adjoint, linearization)
end

function combine_operators(op1, op2)
    @assert op1.target == op2.domain
    operators = vcat(operators(op1), operators(op2))
    return CombinedOperator(op2.domain, op1.target, operators)
end


#Sugar: exp(domain) = PointwiseOperator(domain, exp, exp, exp)

#Operator Interface
#   domain(op) = op.domain
#
#   target(op)
#
#   operators(op)   -> Pointwise: op
#                   -> Linear: op
#                   -> Combined: op.operators
#
#   reverse(op)     -> Pointwise: reverse(ops)
#                   -> Linear: reverse(adjoint(ops))
#
#   getjacobian(op) -> Pointwise: return y -> op'(x) .* y, op(x)
#                   -> Linear: return op(x), op(x)
#


function combine_operators(outer, inner)
    @assert domain(outer) == target(inner)
    operations = _combine_operators(outer, inner)
    return make_Operator(domain(inner), target(outer), operations)
end


struct Model
    domain
    operators
    gradients
    reverse
end

function apply(model::Model, val)
end
function apply(model::Model, val, jacobian)
end

























exp_operator = Operator(domain, exp, exp)
sin(domain) = Operator(domain, sin, cos)
cos(domain) = Operator(domain, cos, -sin)

function getjacobian(val, jacobian, ops, grads)
    gradient = ones(size(val))
    new_val = copy(val)
    for (op, grad) in zip(ops, grads)
        gradient .*= grad.(new_val)
        new_val .= op.(new_val)
    end
    new_jacobian(x) = gradient .* jacobian(x) 
    return new_val, new_jacobian
end


function getjacobian(val, jacobian, ops)
    gradient = ones(size(val))
    new_val = copy(val)
    apply(x) = begin
        new_x = copy(x)
        for op in ops
            new_x .= op(x)
        end
        return new_x
    end
    new_jacobian(x) = apply(jacobian(x))
    return apply(val), new_jacobian
end


