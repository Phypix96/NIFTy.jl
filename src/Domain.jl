#import Base.eltype, Base.size

abstract type AbstractDomain{N} end
abstract type Domain{N}  <: AbstractDomain{N} end
abstract type UnstructuredDomain{N} <: Domain{N} end
abstract type StructuredDomain{N} <: Domain{N} end

dims(::AbstractDomain{N}) where {N} = N
################################################################################
################################################################################
#RGDomain

struct RGDomain{N, H, L} <: StructuredDomain{N}
    _shape::NTuple{N, UInt}
    _lengths::NTuple{N, T} where T
    function RGDomain(shape...; T = Float64, lengths = map(x -> T(1), shape), harmonic = false)
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

struct PowerDomain{N, Dom} <: StructuredDomain{N} where Dom <: RGDomain{_N, true, _L} where {_N, _L}
    _pindices::Array{UInt32, N}
    _kvec::Array{T, 1} where T <: Real
    PowerDomain(domain::RGDomain{N, true, L}, dtype = Float64) where {N, L} = begin
        kvec, karr = get_k_vals(size(domain), distances(domain), dtype)
        pindices = get_pindices(size(domain), kvec, karr)
        typeof(domain) == DataType ? new{N, domain}(pindices, kvec) : new{N, typeof(domain)}(pindices, kvec)
    end
    PowerDomain(shape...; distances = 1 ./ shape, dtype = Float64) = begin
        @assert length(shape) == length(distances)
        kvec, karr = get_k_vals(shape, distances, dtype)
        pindices = get_pindices(shape, kvec, karr)
        domain = RGDomain(shape; distances = distances, harmonic = true, dtype = dtype)
        new{dims(domain), domain}(pindices, kvec)
    end
end


#Get all unique k-values and an array with the corresponding k-values
function get_k_vals(shape, distances, dtype)
    C = CartesianIndices(div.(shape, 2) .+1)
    k_vec = Array{dtype}(undef, length(C))
    k_arr = similar(C, dtype)
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
