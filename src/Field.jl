import Base.eachindex, Base.getindex, Base.randn

include("Domain.jl")

abstract type AbstractField end
abstract type Linearization end


mutable struct Field{domain} <: AbstractField
    val::Array
end

field(domain, val) = begin
    println(shape(domain))
    println(size(val))
    @assert shape(domain) == size(val)
    Field{domain}(val)
end

mutable struct MultiField <: AbstractField
    domains::Tuple
    vals::Tuple
end

#getindex(f::Field, i) = getindex(f.val, i)
getindex(f::Field, i...) = getindex(f.val, i...)
eachindex(f::Field) = eachindex(f.val)
#Use Metaprogramming here!
randn(d::AbstractDomain) = field(d, randn(shape(d)))
rand(d::AbstractDomain) = field(d, rand(shape(d)))
