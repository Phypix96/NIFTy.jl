import Base.size, Base.getindex, Base.setindex!, Base.eltype
import Base.randn, Base.eachindex 

include("Domain.jl")

abstract type AbstractField{T, N} <: AbstractArray{T, N} end
abstract type Linearization end


mutable struct Field{Dom} <: AbstractField{T, N} where {T, N}
    domain::AbstractDomain
    val::Array
end

field(domain, val) = begin
    @assert shape(domain) == size(val)
    Field{typeof(domain)}(domain, val)
    #typeof(domain) == DataType ? Field{domain}(val) : Field{typeof(domain)}(val)
end

mutable struct MultiField <: AbstractField
    domains::Tuple
    vals::Tuple
end

getdomain(f::Field{Dom}) where Dom = Dom

size(f::Field) = size(f.val)
getindex(f::Field, i::Int) = getindex(f.val, i)
getindex(f::Field, I::Vararg{Int, N}) = getindex(f.val, I)
setindex!(f::Field, v, i::Int) = setindex!(f.val, v, i)
setindex!(f::Field, v, I::Vararg{Int, N}) = setindex!(f.val, v, I)
eltype(f::Field) = eltype(f.val)

eachindex(f::Field) = eachindex(f.val)
#Use Metaprogramming here!
randn(d::AbstractDomain) = field(d, randn(shape(d)))
rand(d::AbstractDomain) = field(d, rand(shape(d)))

function get_pindices(ps::Field{PD}) where PD <: PowerDomain
    return ps.domain._pindices
end
