import Base.size, Base.getindex, Base.setindex!, Base.eltype
import Base.randn, Base.eachindex, Base.similar

include("Domain.jl")

abstract type AbstractField{T, N} <: AbstractArray{T, N} end
abstract type AbstractMultiField end
abstract type Linearization end


struct Field{T, N, S, Dom} <: AbstractField{T, N} where Dom <: AbstractDomain
    domain::Dom
    val::Array{T, N}
end

field(domain, val) = begin
    @assert size(domain) == size(val)
    Field{eltype(val), dims(domain), size(domain), typeof(domain)}(domain, val)
end

struct MultiField <: AbstractMultiField
    domains::Tuple
    vals::Tuple
end

getdomain(f::Field{Dom}) where Dom = Dom

#Use Metaprogramming here!
size(f::Field) = size(f.val)
getindex(f::Field, i::Int) = getindex(f.val, i)
getindex(f::Field, I::Vararg{Int, N}) where N = getindex(f.val, I...)
setindex!(f::Field, v, i::Int) = setindex!(f.val, v, i)
setindex!(f::Field, v, I::Vararg{Int, N}) where N = setindex!(f.val, v, I...)
eltype(f::Field) = eltype(f.val)

similar(f::Field) = field(f.domain, similar(f.val))

eachindex(f::Field) = eachindex(f.val)
randn(d::AbstractDomain) = field(d, randn(size(d)...))
rand(d::AbstractDomain) = field(d, rand(size(d)...))

function pindices(ps::Field{T, N, S, PD}) where {T, N, S, PD <: PowerDomain}
    return ps.domain._pindices
end
