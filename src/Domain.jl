import Base.eltype, Base.size

abstract type AbstractDomain end
abstract type UnstructuredDomain <: AbstractDomain end
abstract type StructuredDomain <: AbstractDomain end

struct RGDomain{N, S, D, H, T} <: StructuredDomain
    #Those fields only exist to enforce the length and data type of shape and distances
    shape::NTuple{N,UInt32}
    distances::NTuple{N,T}
    RGDomain(shape...; distances = 1 ./ shape, harmonic = false, dtype = Float64) = new{length(shape), shape, distances, harmonic, dtype}(shape, distances)
    #RGDomain(shape...; distances = 1., harmonic = true, dtype = Float64) = new{length(shape), shape, distances, harmonic, dtype}(shape, distances)
end

shape(::RGDomain{N, S, D, H, T}) where {N, S, D, H, T} = S
distances(::RGDomain{N, S, D, H, T}) where {N, S, D, H, T} = D
is_harmonic(::RGDomain{N, S, D, H, T}) where {N, S, D, H, T} = H
eltype(::RGDomain{N, S, D, H, T}) where {N, S, D, H, T} = T
is_complex(::RGDomain{N, S, D, H, T}) where {N, S, D, H, T} = T <: Complex

function get_codomain(domain::RGDomain) 
    return RGDomain(shape(domain)...;
                    distances = shape(domain) .* distances(domain),
                    harmonic = !is_harmonic(domain),
                    dtype = eltype(domain))
end

################################################################################
################################################################################
#Power Domain

struct PowerDomain{Dom} <: StructuredDomain where Dom <: RGDomain
    #_pindices::Vector{length(shape(Dom)), UInt32}
    #_kvec::Vector{1, eltype(Dom)}
    _pindices::Array{N, UInt32} where N
    _kvec::Array{1, T} where T
    PowerDomain(domain::RGDomain{N, S, D, true, T}) where {N, S, D, T} = begin
        kvec, karr = get_k_vals(shape(domain), distances(domain), eltype(domain))
        pindices = get_pindices(shape(domain), kvec, karr)
        typeof(domain) == DataType ? new{domain}(pindices, kvec) : new{typeof(domain)}(pindices, kvec)
    end
end

#PowerDomain(shape...; distances = 1 ./ shape, dtype = Float64) = begin
#    @assert length(shape) == length(distances)
#    kvec, karr = get_k_vals(shape, distances, dtype)
#    pindices = get_pindices(shape, kvec, karr)
#    domain = RGDomain(shape; distances = distances, harmonic = true, dtype = dtype)
#    new{domain}(pindices, kvec)
#end

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


shape(dom::PowerDomain) = size(dom._kvec)
