using NIFTy
using Plots

dom = RGDomain(200,200, harmonic = true)
pdom = PowerDomain(dom)

f1 = randn(size(dom))
f2 = randn(size(dom))
kernel1= map(x -> 1/(1+x)^1.4, pdom._kvec)
kernel2 = map(x -> 1/(1+x)^3, pdom._kvec)

distributor = DistributionOperator((dom, pdom), dom, pdom._pindices)
ht = hartley(dom)
op = exp(ht(distributor))

res1 = op(f1, kernel1)
res2 = op(f2, kernel2)

heatmap(res1 + 50*res2)
