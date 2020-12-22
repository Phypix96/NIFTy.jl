#import Base: size, getindex, setindex!, IndexStyle, randn, similar

#using .Domain
#include("Domain.jl")

################################################################################
################################################################################
#Field type
abstract type AbstractField{T, N, Dom} <: AbstractArray{T, N} where Dom <: AbstractDomain end
struct Field{T, N, Dom} <: AbstractField{T, N, Dom}
#struct Field{T, N, Dom} <: AbstractArray{T, N} where Dom <: AbstractDomain
    domain::Dom
    val::AbstractArray{T, N}
    #val::A{T, N} where A <: AbstractArray{T, N}
    Field(domain, val) = begin
        @assert size(domain) == size(val)
        new{eltype(val), dims(domain), typeof(domain)}(domain, val)
    end
end

#Implement correct broadcasting behaviour
#Pointwise operations on a field should always default to preserving the domain
#Binary broadcasting is only well defined for fields with matching domains
struct FieldStyle{N, Dom} <: Broadcast.AbstractArrayStyle{N} end
Base.BroadcastStyle(::Type{<:Field{T, N, Dom}}) where {T, N, Dom} = FieldStyle{N, Dom}()
Base.BroadcastStyle(::FieldStyle{N, Dom}, ::FieldStyle{N, Dom}) where {N, Dom} = FieldStyle{N, Dom}()
Base.BroadcastStyle(::FieldStyle, ::FieldStyle) = error("Domain mismatch")
#Base.BroadcastStyle(::Broadcast.Style{<:Number}, ::FieldStyle{N, Dom}) where {N, Dom} = FieldStyle{N, Dom}()
#Base.BroadcastStyle(::FieldStyle{N, Dom}, ::Broadcast.Style{<:Number}) where {N, Dom} = FieldStyle{N, Dom}()

similar(f::Field) = similar(f, eltype(f), size(f))
similar(f::Field, ::Type{T}, ::Tuple{Vararg{Int64, N}}) where {T, N} = Field(f.domain, similar(f.val, T))
function similar(bc::Broadcast.Broadcasted{FieldStyle{N, Dom}}, ::Type{ElType}) where {N, Dom, ElType}
    domain = getdomain(bc)
    valtype = getvaltype(bc)
    return Field(domain, similar(valtype, axes(bc)))
end

getdomain(bc::Base.Broadcast.Broadcasted) = getdomain(bc.args[1])
getdomain(f::Field) = f.domain
getvaltype(bc::Base.Broadcast.Broadcasted) = getvaltype(bc.args[1])
getvaltype(f::Field) = typeof(f.val)

shape(f::Field) = shape(f.domain)
size(f::Field) = Int.(size(f.val))
IndexStyle(::Type{<:Field}) = IndexLinear()

getindex(f::Field, i::Integer) = getindex(f.val, i)
getindex(f::Field, I::Vararg{Integer, N}) where N = getindex(f.val, I...)
setindex!(f::Field, v, i::Integer) = setindex!(f.val, v, i)
setindex!(f::Field, v, I::Vararg{Integer, N}) where N = setindex!(f.val, v, I...)


randn(d::AbstractDomain) = Field(d, randn(size(d)...))
randn(::Type{T}, d::AbstractDomain) where T = Field(d, randn(T, size(d)...))
randn(::Type{T}, d::AbstractDomain) where T <: Integer = Field(d, round.(T, randn(size(d)...)))


pindices(ps::Field{T, N, PD}) where {T, N, PD <: PowerDomain} = ps.domain._pindices

#MultiField: named tuples of fields



struct Linearization{T, N, Dom} <: AbstractField{T, N, Dom}
    val::Field{T, N, Dom}
    jacobian
end
#
#struct LinearizationStyle{N, Dom} <: Broadcast.AbstractArrayStyle{N} end
#Base.BroadcastStyle(::Type{<:Linearization{T, N, Dom}}) where {T, N, Dom} = LinearizationStyle{N, Dom}()
#Base.BroadcastStyle(::LinearizationStyle{N, Dom}, ::LinearizationStyle{N, Dom}) where {N, Dom} = LinearizationStyle{N, Dom}()
#Base.BroadcastStyle(::LinearizationStyle, ::LinearizationStyle) = error("Domain mismatch")
#
#function similar(bc::Broadcast.Broadcasted{LinearizationStyle{N, Dom}}, ::Type{ElType}) where {N, Dom, ElType}
#    val = getval(bc)
#    jacobian = getjacobian(bc)
#    new_jacobian(x) = val .* jacobian(x)
#    return Linearization(similar(val), new_jacobian)
#end
