using Plots
using Random, Distributions
using LinearAlgebra
using WAV
using DSP

purples = ["#0a0612", "#392249", "#482980", "#673ab7",
    "#7a52aa", "#9779bd", "#b59fd0", "#d3c5e3"]

function getharmonics(f1, n)
    harmonics = zeros(n)
    for i = 1:length(harmonics)
        harmonics[i] = i * f1
    end
    return harmonics
end

function average(x)
    return sum(x) / length(x)
end

function harmonic(t; f1=1, alist=1, Φlist=0)
    f = zeros(length(alist))
    for i = 1:length(alist)
        f[i] = i * f1
    end
    value = alist .* cos.(2 .* pi .* (f .* t .- Φlist))
    return sum(value)
end

function cosines(t; flist=1, alist=1, Φlist=0)
    value = alist .* cos.(2 .* pi .* (flist .* t .- Φlist))
    return sum(value)
end

function writetowav(f; duration=3.0, filename="1b.wav", createfile=true, playfile=true)
    # This function writes the function f to a wav file, plays it, and returns Vectors for plotting purposes
    fs = 44100
    t = 0.0:1/fs:prevfloat(duration)
    signal = 0.1 * f.(t)
    if createfile
        wavwrite(signal, string("output_sounds/", filename), Fs=fs)
    end
    if playfile
        wavplay(signal, fs)
    end
    return (signal, t) # for plotting purposes
end

function signalenergy(signal)
    return norm(signal, 2)^2
end

function signalpower(signal)
    return signalenergy(signal) / length(signal)
end

function createnoise(; noisetype="Gaussian", duration=3.0, filename="noise.wav", createfile=true, playfile=true)
    # This function writes noise to a wav file, plays it, and returns Vectors for plotting purposes
    fs = 44100
    clipped = false
    t = 0.0:1/fs:prevfloat(duration)
    signal = []
    if noisetype == "Gaussian"
        signal = 0.1 * randn(length(t))
    elseif noisetype == "Laplacian"
        signal = 0.709 * 0.1 * rand(Laplace(), length(t))
    elseif noisetype == "Uniform"
        signal = 1.733 * 0.1 * rand(-1.0:eps():1.0, length(t))
    else
        throw(ArgumentError("Try noisetype=Gaussian, Laplacian, or Uniform"))
    end
    println("Signal Power: ", signalpower(signal))
    # avoiding clipping:
    for i in 1:length(t)
        if signal[i] >= 1.0
            signal[i] = prevfloat(1.0)
            clipped = true
        end
    end
    if createfile
        wavwrite(signal, string("output_sounds/", filename), Fs=fs)
    end
    if playfile
        wavplay(signal, fs)
    end
    println("Distribution clipped? ", clipped)
    return (signal, t) # for plotting purposes
end

function gabor(t; a=1.0, σ=1, f=1.0, Φ=0.0)
    # Note: the Φ here is different than the Φ in harmonic() and cosine().  It has an effective range of [0, 2π)
    return a * exp((-t^2) / (2 * (σ^2))) * cos(2 * pi * f * t + Φ)
end

function filterIIR(x, a, b)
    y = copy(x)
    for i in 1:length(y)
        for j in 1:length(b)
            y[i] += b[j] * x[i-j+1]
        end
        for k in 1:length(a)
            if (i - k) >= 1
                y[i] -= a[k] * y[i-k]
            end
        end
    end
    return y
end

function createfilterednoise(f, σ; noisetype="Gaussian", filtertype="IIR", duration=3.0, filename="filtered_noise.wav", createfile=true, playfile=true)
    fs = 44100
    signal, t = createnoise(noisetype="Gaussian", filename="doesntmatter.wav", duration=duration, createfile=false, playfile=false)
    filteredsignal = []
    if filtertype == "IIR"
        filteredsignal = 0.1 * filterIIR(signal, [-1.702, 0.81], [0.063])
    else
        ArgumentError("other filtering remains not implemented; try filtertype = IIR ")
        #=envelopet = [reverse(t); t[2:end]] 
        gaborenvelope = gabor.(envelopet; f=f, σ=σ)
        filterbleed = Int64(ceil((length(gaborenvelope)/2)-1))
        filteredsignal = 0.1 * conv(signal, gaborenvelope)=#
    end

    # avoiding clipping:
    for i in 1:length(t)
        if filteredsignal[i] >= 1.0
            filteredsignal[i] = prevfloat(1.0)
        end
    end

    if createfile
        wavwrite(signal, string("output_sounds/", filename), Fs=fs)
    end
    if playfile
        # wavplay(filteredsignal[begin+filterbleed:end-filterbleed], fs)
        wavplay(filteredsignal, fs)
    end
    return (filteredsignal, t)
end

function convolve(x, y)
    # y is 'scanning' x in a sense
    x_new = [zeros(length(y)); copy(x); zeros(length(y))]
    result = zeros(length(x) + length(y))
    for n = 1:length(result)
        for k = 1:length(y)
            result[n] += y[k] * x_new[length(y)+n-k+1]
        end
    end
    return result[begin:end-1]
end

function crosscorr(g, x; normalize=true)
    # g is 'scanning' x in a sense
    x_new = [zeros(length(g)); copy(x); zeros(length(g))]
    result = zeros(length(g) + length(x))
    for i = 1:length(result)
        for j = 1:length(g)
            result[i] += x_new[i+j-1] * g[j]
        end
    end
    if normalize
        return result[begin+1:end] / (norm(x) * norm(g))
    else
        return result[begin+1:end]
    end
end

function autocorr(x; normalize=true)
    crosscorr(x, x, normalize=normalize)
end

function boxwaveform(signal, t; n=500)
    new_signal = []
    new_t = []
    gap = length(signal) / n
    indices = Int64.(round.(1:gap:length(signal)))
    for i in 1:(length(indices)-1)
        push!(new_signal, average(signal[indices[i]:indices[i+1]])) # sum(signal[i:i+Int64(round(gap))-1])
        push!(new_t, t[indices[i]])
    end
    return (new_signal, new_t)
end

function boxwaveform(signal; n=500)
    new_signal = []
    gap = length(signal) / n
    indices = Int64.(round.(1:gap:length(signal)))
    for i in 1:(length(indices)-1)
        push!(new_signal, average(signal[indices[i]:indices[i+1]])) # sum(signal[i:i+Int64(round(gap))-1])
    end
    return new_signal
end

function estimatepitch(waveform; thresh=0.0500)
    fs = 44100

    crossings = []
    peakcenters = []
    withinpeak = false
    for i in 2:length(waveform)
        if (waveform[i-1] <= thresh && waveform[i] > thresh)
            withinpeak = true
            push!(crossings, i)
        end
        if (waveform[i-1] > thresh && waveform[i] <= thresh)
            withinpeak = false
            push!(crossings, i)
            push!(peakcenters, crossings[end-1] + ((crossings[end] + crossings[end-1]) / 2))
        end
    end
    return 44100 / (sum(abs.(peakcenters[begin:end-1] .- peakcenters[begin+1:end])) / (length(peakcenters) - 1) / 2)
end;

function estimatedelay(signal1, signal2; fs=44100)
    correlation = crosscorr(signal1, signal2, normalize=true)

    max = -1.0
    maxindex = 1  # this will be the index of the signal with the highest correlation
    for i in 1:length(correlation)
        if correlation[i] >= max
            max = correlation[i]
            maxindex = i
        end
    end
    delayinsamples = maxindex - length(signal1)
    delayinseconds = delayinsamples / fs
    return delayinseconds
end