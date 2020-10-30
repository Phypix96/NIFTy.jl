using Plots, FFTW

include("Domain.jl")
include("Field.jl")
include("Distributor.jl")

dom = RGDomain(100,100)

hdom = get_codomain(dom)
pdom = PowerDomain(hdom)

f = randn(hdom)
fcorr = field(hdom, similar(f.val))
kernel = field(pdom, 1 ./ (1 .+ pdom._kvec).^1.5)

PowerDistribute!(f, kernel, fcorr)

#heatmap(real(fft(fcorr.val)))

