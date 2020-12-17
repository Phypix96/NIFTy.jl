using NIFTy
using Plots

dom = RGDomain(200,200, harmonic = true)
pdom = PowerDomain(size(dom)...)

f1 = randn(size(dom))
f2 = randn(size(dom))
kernel1= map(x -> 1/(1+x)^1.5, pdom._kvec)
kernel2 = map(x -> 1/(1+x)^3, pdom._kvec)

distributor = DistributionOperator(dom, dom, pdom._pindices)
ht = hartley(dom)
op = exp(getcodomain(dom))

res1 = distributor(f1, kernel1)
res2 = distributor(f2, kernel2)
#apply!(ht, res)
res1 = op(ht(res1))
res2 = op(ht(res2))
heatmap(res1 + 100*res2)
