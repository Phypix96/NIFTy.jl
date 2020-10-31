import Base.size, Base.length, Base.getindex, Base.setindex!, Base.eltype, Base.iterate
import Base.randn, Base.eachindex, Base.similar

include("Domain.jl")

#abstract type AbstractField{T, N} <: AbstractArray{T, N} end
abstract type AbstractField{T, N} end
abstract type AbstractMultiField end
abstract type Linearization end


struct Field{T, N, Dom} <: AbstractField{T, N} where Dom <: AbstractDomain
    domain::Dom
    val::Array{T, N}
end

field(domain, val) = begin
    @assert size(domain) == size(val)
    Field{eltype(val), dims(domain), typeof(domain)}(domain, val)
end

struct MultiField <: AbstractMultiField
    domains::Tuple
    vals::Tuple
end

getdomain(f::Field) = f.domain


#Use Metaprogramming here!
size(f::Field) = size(f.val)
length(f::Field) = prod(size(f))
getindex(f::Field, i::Integer) = getindex(f.val, i)
getindex(f::Field, I::Vararg{Integer, N}) where N = getindex(f.val, I...)
setindex!(f::Field, v, i::Integer) = setindex!(f.val, v, i)
setindex!(f::Field, v, I::Vararg{Integer, N}) where N = setindex!(f.val, v, I...)
eltype(f::Field) = eltype(f.val)

similar(f::Field) = field(f.domain, similar(f.val))
similar(f::Field, dom::AbstractDomain) = field(dom, similar(f.val))

eachindex(f::Field) = eachindex(f.val)
randn(d::AbstractDomain) = field(d, randn(size(d)...))
rand(d::AbstractDomain) = field(d, rand(size(d)...))

function pindices(ps::Field{T, N, PD}) where {T, N, PD <: PowerDomain}
    return ps.domain._pindices
end
