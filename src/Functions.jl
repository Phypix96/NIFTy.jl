include("Domain.jl")
include("Field.jl")

################################################################################
################################################################################
#Fourier and Hartley transformation
function fft(f::Field{T, N, Dom}) where {T, N, Dom <: RGDomain}
    dom = getdomain(f)
    if !isharmonic(dom)
        error("Fourier Transform requires harmonic domain")
    end
    new_dom = getcodomain(dom)
    val = fft(f.val)
    return field(new_dom, val)
end

function ifft(f::Field{T, N, Dom}) where {T, N, Dom <: RGDomain}
    dom = getdomain(f)
    if isharmonic(dom)
        error("Inverse Fourier Transform requires real domain")
    end
    new_dom = getcodomain(dom)
    val = ifft(f.val)
    return field(new_dom, val)
end

function hartley(f::Array)
    return FFTW.r2r(f, FFTW.DHT)
end

function ihartley(f::Array)
    return hartley(f) ./ size(f)
end

