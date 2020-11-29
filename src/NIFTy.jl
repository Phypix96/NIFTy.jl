module NIFTy

import Base: size, length, eltype, getindex, setindex!, IndexStyle, randn, similar

import FFTW: fft, ifft

include("Domain.jl")
include("Field.jl")
include("Functions.jl")
include("Distributor.jl")

export AbstractDomain, Domain, UnstructuredDomain, StructuredDomain,
       RGDomain, PowerDomain, DomainTuple,
       rgdomain,
       lengths, isharmonic, distances, dims, shape,
       getcodomain, 
       domains,
       #Field.jl
       Field,
       similar, getdomain,
       #Functions.jl
       hartley, ihartley


end
