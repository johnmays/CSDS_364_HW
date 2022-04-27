function plotfourierbasis(;k, N)
    # plotting the unit circle
    xₜ(t) = sin(t)
    yₜ(t) = cos(t)
    plot(xₜ, yₜ, 0, 2π, leg=false, color = "black", 
    size = (400,400),framestyle = :origin,
    title="Fourier Basis", xlabel = "Re", ylabel = "Im",
    titlefontsize = 12, xguidefontsize = 8, yguidefontsize = 8, tickfontsize = 6,
    xlim = (-1.5, 1.5), ylim = (-1.5, 1.5)
    )
    # plotting the vector
    ω = (2*pi*k)/N
    basis = exp(im*ω)
    plot!([0,real(basis)], [0,imag(basis)], color = purples[5], lw=2)
    # print(round.(basis, digits =3)) 
    plot!()
end

function plotfourierbasis(N; k=1/N)
    # plotting the unit circle
    xₜ(t) = sin(t)
    yₜ(t) = cos(t)
    plot(xₜ, yₜ, 0, 2π, leg=false, color = "black", 
    size = (400,400),framestyle = :origin,
    title="Fourier Basis", xlabel = "Re", ylabel = "Im",
    titlefontsize = 12, xguidefontsize = 8, yguidefontsize = 8, tickfontsize = 6,
    xlim = (-1.5, 1.5), ylim = (-1.5, 1.5)
    )
    # plotting the vectors
    basis = zeros(ComplexF64,N)
    for k = 0:(N-1)
        ω = (2*pi*k)/N
        basis[i] = exp(im*ω)
        plot!([0,real(basis[i])], [0,imag(basis[i])], color = purples[5], lw=3)
          
    end
    # print(round.(basis, digits =3)) 
    plot!()
end