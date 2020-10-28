include("Domain.jl")
include("Field.jl")
include("Distributor.jl")

dom = RGDomain(100,100)
hdom = get_codomain(dom)
pdom = PowerDomain(hdom)

f = randn(hdom)
fcorr = field(dom, similar(f.val))
kernel = field(pdom, 1 ./ pdom._kvec.^2)

PowerDistribute!(f, kernel, fcorr)



