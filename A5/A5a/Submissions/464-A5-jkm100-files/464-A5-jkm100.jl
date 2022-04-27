using Plots
using Random, Distributions
using LinearAlgebra
using DSP
using FFTW
using Images

purples = ["#0a0612", "#392249", "#482980", "#673ab7",
    "#7a52aa", "#9779bd", "#b59fd0", "#d3c5e3"]

function plotfourierbasis(N; k=nothing)
    # plotting the unit circle
    xₜ(t) = sin(t)
    yₜ(t) = cos(t)
    plot(xₜ, yₜ, 0, 2π, leg=false, color="black",
        size=(400, 400), framestyle=:origin,
        title="Fourier Basis", xlabel="Re", ylabel="Im",
        titlefontsize=12, xguidefontsize=8, yguidefontsize=8, tickfontsize=6,
        xlim=(-1.5, 1.5), ylim=(-1.5, 1.5)
    )
    # plotting the vectors
    if k isa Number
        ω = (2 * pi * k) / N
        basis = exp(im * ω)
        plot!([0, real(basis)], [0, imag(basis)], color=purples[5], lw=2)
    else
        basis = zeros(ComplexF64, N)
        for k = 0:(N-1)
            ω = (2 * pi * k) / N
            basis[k+1] = exp(im * ω)
            plot!([0, real(basis[k+1])], [0, imag(basis[k+1])], color=purples[5], lw=1)

        end
    end
    # print(round.(basis, digits =3)) 
    plot!()
end

function w(n, k, N)
    ω = (2 * pi * k) / N
    return exp(im * ω * n)
end

function plotw(k, N)
    basis = zeros(ComplexF64, N)
    for n = 0:(N-1)
        basis[n+1] = w(n, k, N)
    end
    real_plot = plot(real(basis), line=:stem, ylabel="Re", label=false,
        color=purples[3])
    imag_plot = plot(imag(basis), line=:stem, ylabel="Im", label=false,
        color=purples[3])

    plot(real_plot, imag_plot,
        layout=(2, 1), size=(500, 250), xguidefontsize=8, yguidefontsize=8
    )
end

function fourier_matrix(N::Int64)
    A = zeros(ComplexF64, (N, N))
    for row = 1:N
        n = row - 1
        for col = 1:N
            k = col - 1
            A[row, col] = w(n, k, N)
        end
    end
    return A
end

function myft(y)
    N = length(y)
    return fourier_matrix(N) * y
end

function w2D(k, j, n, m, N)
    return cos((2 * pi * k * n) / N) + cos((2 * pi * j * m) / N)
end

function make2Dcell(k, j, N)
    cell = zeros(Float64, (N, N))
    for row in 1:N
        n = row - 1
        for col in 1:N
            m = col - 1
            cell[row, col] = w2D(k, j, n, m, N)
        end
    end
    return cell
end