#import Base: size, getindex, setindex!, IndexStyle, randn, similar

#using .Domain
#include("Domain.jl")

################################################################################
################################################################################
#Field type
struct Field{T, N, Dom} <: AbstractArray{T, N} where Dom <: AbstractDomain
    domain::Dom
    val::Array{T, N}
    Field(domain, val) = begin
        @assert size(domain) == size(val)
        new{eltype(val), dims(domain), typeof(domain)}(domain, val)
    end
end

struct MultiField
    domain::Union{Tuple, NamedTuple}
    val::Union{Tuple, NamedTuple}
end

function MultiField(fs::Union{Field, MultiField}...)
    domains = Tuple(f.domain for f in fs)
    vals = Tuple(f.val for f in fs)
    return MultiField(domains, vals)
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
    return Field(domain, similar(Array{ElType}, axes(bc)))
end

getdomain(bc::Base.Broadcast.Broadcasted) = getdomain(bc.args...)
#TODO assure, that all args have the same domain?
getdomain(args::Vararg{Field}) = getdomain(args[1])
getdomain(f::Field) = f.domain

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


