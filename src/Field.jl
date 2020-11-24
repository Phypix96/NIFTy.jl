import Base.size, Base.getindex, Base.setindex!, Base.IndexStyle
import Base.randn, Base.similar

include("Domain.jl")

struct Field{A, T, N, Dom} <: AbstractArray{T, N} where Dom <: AbstractDomain
    domain::Dom
    val::Array{T, N}
end

field(domain, val, adjoint=false) = begin
    @assert size(domain) == size(val)
    Field{adjoint, eltype(val), dims(domain), typeof(domain)}(domain, val)
end

struct FieldStyle{A, N, Dom} <: Broadcast.AbstractArrayStyle{N} end
Base.BroadcastStyle(::Type{<:Field{A, T, N, Dom}}) where {A, T, N, Dom} = FieldStyle{A, N, Dom}()
Base.BroadcastStyle(::FieldStyle{A, N, Dom}, ::FieldStyle{A, N, Dom}) where {A, N, Dom} = FieldStyle{A, N, Dom}()
#Base.BroadcastStyle(::Broadcast.Style{<:Number}, ::FieldStyle{A, N, Dom}) where {A, N, Dom} = FieldStyle{A, N, Dom}()
#Base.BroadcastStyle(::FieldStyle{A, N, Dom}, ::Broadcast.Style{<:Number}) where {A, N, Dom} = FieldStyle{A, N, Dom}()

function similar(bc::Broadcast.Broadcasted{FieldStyle{A, N, Dom}}, ::Type{ElType}) where {A, N, Dom, ElType}
    domain = getdomain(bc)
    return field(domain, similar(Array{ElType}, axes(bc)), A)
end

getdomain(bc::Base.Broadcast.Broadcasted) = getdomain(bc.args)
getdomain(args::Tuple) = getdomain(args[1])
getdomain(f::Field) = f.domain

is_adjoint(f::Field{A}) where A = A

#Use Metaprogramming here!
size(f::Field) = Int.(size(f.val))
IndexStyle(::Type{<:Field}) = IndexLinear()

getindex(f::Field, i::Integer) = getindex(f.val, i)
getindex(f::Field, I::Vararg{Integer, N}) where N = getindex(f.val, I...)
setindex!(f::Field, v, i::Integer) = setindex!(f.val, v, i)
setindex!(f::Field, v, I::Vararg{Integer, N}) where N = setindex!(f.val, v, I...)


similar(f::Field) = similar(f, eltype(f), size(f))
similar(f::Field, ::Type{T}, ::Tuple{Vararg{Int64, N}}) where {T, N} = field(f.domain, similar(f.val, T))
#similar(f::Field) = field(f.domain, similar(f.val))
#similar(f::Field, dom::AbstractDomain) = field(dom, similar(f.val))


randn(d::AbstractDomain) = field(d, randn(size(d)...))
rand(d::AbstractDomain) = field(d, rand(size(d)...))

function pindices(ps::Field{T, N, PD}) where {T, N, PD <: PowerDomain}
    return ps.domain._pindices
end

#MultiField: named tuples of fields
