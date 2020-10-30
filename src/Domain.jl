import Base.eltype, Base.size

abstract type AbstractDomain end
abstract type UnstructuredDomain <: AbstractDomain end
abstract type StructuredDomain <: AbstractDomain end

struct RGDomain{N, H, L} <: StructuredDomain
    _shape::NTuple{N, UInt}
    _lengths::NTuple{N, T} where T
    RGDomain(shape...; lengths = map(x -> 1., shape), harmonic = false) = begin
        new{length(shape), harmonic, lengths}(shape, lengths)
    end
end

#function RGDomain(shape...; distances = 1 ./ shape, harmonic = false)
    #RGDomain(shape; lengths = shape .* distances, harmonic)
#end
size(domain::RGDomain) = domain._shape
lengths(::RGDomain{N, H, L}) where {N, H, L} = L
dims(::RGDomain{N, H, L}) where {N, H, L} = N
is_harmonic(::RGDomain{N, H, L}) where {N, H, L} = H
distances(domain::RGDomain) = lengths(domain) ./ size(domain)

function get_codomain(domain::RGDomain) 
    return RGDomain(size(domain)...; lengths = size(domain) ./ lengths(domain), harmonic = !is_harmonic(domain))
end

################################################################################
################################################################################
#Power Domain

struct PowerDomain{N, Dom} <: StructuredDomain where Dom <: RGDomain{N, true, L} where {N <: Integer, L <: DataType}
    _pindices::Array{UInt32, N}
    _kvec::Array{T, 1} where T <: Real
    PowerDomain(domain::RGDomain{N, true, L}, dtype = Float64) where {N, L} = begin
        kvec, karr = get_k_vals(size(domain), distances(domain), dtype)
        pindices = get_pindices(size(domain), kvec, karr)
        typeof(domain) == DataType ? new{N, domain}(pindices, kvec) : new{N, typeof(domain)}(pindices, kvec)
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


size(dom::PowerDomain) = size(dom._kvec)
dims(::PowerDomain) = 1
