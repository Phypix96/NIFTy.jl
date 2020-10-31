#using Plots
using FFTW

include("Domain.jl")
include("Field.jl")
include("Distributor.jl")

dom = RGDomain(512,512)

hdom = get_codomain(dom)
pdom = PowerDomain(hdom)

f = randn(hdom)
fcorr = similar(f, hdom)
kernel = field(pdom, 1 ./ (1 .+ pdom._kvec).^1.5)

fcorr = PowerDistribute(f, kernel);

corr = real(fft(fcorr))
#heatmap(real(fft(fcorr.val)))

