module NIFTy

import Base: size, length, eltype, getindex, setindex!, IndexStyle, randn, similar
import Base: adjoint
import Base: exp, log, sin, cos

import FFTW: fft, ifft

using LoopVectorization

include("Domain.jl")
include("Field.jl")
include("Operator.jl")
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
       hartley, ihartley, 
       #Operator.jl
       AbstractOperator, LinearOperator, PointwiseOperator, OperatorChain, Operator,
       domain, target, adjoint,
       combine_operators, apply!, apply,
       add, scale,
       #Distributor.jl
       DistributionOperator, PowerDistributor
end
