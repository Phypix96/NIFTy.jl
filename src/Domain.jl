#import Base.eltype, Base.size

abstract type AbstractDomain{N} end
abstract type Domain{N}  <: AbstractDomain{N} end
abstract type UnstructuredDomain{N} <: Domain{N} end
abstract type StructuredDomain{N, H} <: Domain{N} where H end

dims(::AbstractDomain{N}) where {N} = N
################################################################################
################################################################################
#RGDomain

struct RGDomain{N, H, L} <: StructuredDomain{N, H}
    _shape::NTuple{N, UInt}
    _lengths::NTuple{N, T} where T
    RGDomain(shape...;
             T = Float64,
             harmonic = false,
             lengths = harmonic ? T.(shape) : ones(T, length(shape))) = begin
        lengths = tuple(T.(lengths)...)
        new{length(shape), harmonic, lengths}(shape, lengths)
    end
end

shape(domain::RGDomain) = size(domain)
#size(domain::RGDomain) = domain._shape
size(domain::RGDomain) = Int.(domain._shape)
length(domain::RGDomain) = prod(size(domain))
lengths(::RGDomain{N, H, L}) where {N, H, L} = L
isharmonic(::RGDomain{N, H, L}) where {N, H, L} = H
distances(domain::RGDomain) = lengths(domain) ./ size(domain)

function getcodomain(domain::RGDomain) 
    return RGDomain(size(domain)...; lengths = size(domain) ./ lengths(domain), harmonic = !isharmonic(domain))
end

################################################################################
################################################################################
#Power Domain

struct PowerDomain{N, Dom} <: StructuredDomain{N, true} where Dom <: RGDomain{_N, true, _L} where {_N, _L}
    _codomain::StructuredDomain{N, true}
    _pindices::Array{UInt32, N}
    _kvec::Array{T, 1} where T <: Real
    PowerDomain(domain::RGDomain{N, true, L}, T = Float64) where {N, L} = begin
        kvec, karr = get_k_vals(size(domain), distances(domain), T)
        pindices = get_pindices(size(domain), kvec, karr)
        new{N, typeof(domain)}(domain, pindices, kvec)
    end
end

function PowerDomain(shape::Vararg{Integer}; T = Float64, lengths = T.(shape))
    @assert length(shape) == length(lengths)
    distances = lengths ./ shape
    domain = RGDomain(shape...; T = T, lengths = lengths, harmonic = true)
    return PowerDomain(domain, T)
end


#Get all unique k-values and an array with the corresponding k-values
function get_k_vals(shape, distances, T)
    C = CartesianIndices(div.(shape, 2) .+1)
    k_vec = Array{T}(undef, length(C))
    k_arr = similar(C, T)
    for (i, I) in enumerate(C)
        k = sqrt(sum(((Tuple(I) .- 1) .*distances).^2))
        k_vec[i] = k
        k_arr[I] = k
    end
    unique!(k_vec)
    sort!(k_vec)
    return k_vec, k_arr
end

#Get an array of the same size as the domain with each value being the index of the corresponding value in k-vector
function get_pindices(shape, kvec, karr)
    C = CartesianIndices(shape)
    pindex = similar(C, UInt32)
    p_part = indexin(karr, kvec)
    small_s = Tuple(shape)./2 .+ 1
    for I in C
        pindex[I] = p_part[CartesianIndex(Int.(small_s .- abs.(Tuple(I).-small_s)))]
    end
    return pindex
end

size(dom::PowerDomain) = size(dom._kvec)
dims(::PowerDomain) = 1
codomain(dom::PowerDomain) = dom._codomain

################################################################################
################################################################################
#DomainTuple

struct DomainTuple{N} <: AbstractDomain{N}
    domains::Tuple{Vararg{<:Domain}}
    DomainTuple(domains...) = begin
        N = mapreduce(dims, +, domains)
        new{N}(domains)
    end
end

domains(dom::DomainTuple) = dom.domains
shape(dom::DomainTuple) = map(size, domains(dom))
size(dom::DomainTuple) = tuple(mapreduce(x -> [size(x)...], vcat, domains(dom))...)
length(dom::DomainTuple) = prod(size(dom))
#length(dom::DomainTuple) = prod(length(domains))
