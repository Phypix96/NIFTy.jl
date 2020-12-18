using NIFTy
using Plots

dom = RGDomain(500,500, harmonic = true)
pdom = PowerDomain(dom)

f1 = randn(size(dom))
f2 = randn(size(dom))
kernel1= map(x -> 1/(1+x)^1.4, pdom._kvec)
kernel2 = map(x -> 1/(1+x)^3, pdom._kvec)

distributor = DistributionOperator((dom, pdom), dom, pdom._pindices)
op1 = exp(hartley(distributor))
op2 = scale(target(op1), 80)(op1)
full_op = op1 * op2
heatmap(full_op((f1, kernel1), (f2, kernel2)))

