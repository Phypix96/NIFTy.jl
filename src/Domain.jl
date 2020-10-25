abstract type AbstractDomain end

abstract type UnstructuredDomain <: AbstractDomain end

abstract type StructuredDomain <: AbstractDomain end

struct RGDomain <: StructuredDomain
    shape
    distances
    harmonic::Bool
    dtype::DataType
    RGDomain(shape...; distances = 1 ./ shape, harmonic = false, dtype = Float64) = new(shape, distances, harmonic, dtype)
end

################################################################################
################################################################################
#Power Domain
struct PowerDomain
    shape
    distances
    dtype::DataType
    _pindices
    _kvec
    PowerDomain(shape...; distances = 1 ./ shape, dtype = Float64) = begin
        @assert length(shape) == length(distances)
        kvec, karr = get_k_vals(shape, distances, dtype)
        pindices = get_pindices(shape, kvec, karr)
        new(shape, distances, dtype, pindices, kvec)
    end
end

PowerDomain(domain::RGDomain) = PowerDomain(domain.shape; distances = domain.distances, dtype = domain.dtype)

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



#field[I] = kvec[pindex[I]]
