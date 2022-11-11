using KolmogorovZurbenkoFilters
using Test
using Distributions

@testset "d and dprime (w=$(w))" for w in [3,4,5]
    y_sample = collect(1:20)  # line with unit slope
    d_sample, dprime_sample = KolmogorovZurbenkoFilters.R_differenced(y_sample, w)
    @test d_sample[1:w] == reverse(d_sample[end-w+1:end])
    @test all(d_sample[w+1:end-w] .== w*2)
    @test all(dprime_sample[1:w] .== 1)
    @test all(dprime_sample[end-w:end] .== -1)
    @test all(dprime_sample[w+1:end-w-1] .== 0)
end

@testset "moving average" begin
    N = 40
    N_step = 10
    i = 1:N
    y_ref = -1.0e-1 .* [(i_ + i_ % N_step)/N for i_ in i]
    y_noisy = y_ref .+ 0.05 .* rand(Normal(), N)

    y_mavg = [KolmogorovZurbenkoFilters.mavg1d(y_noisy, i, 9) for i in 1:N]
    # drop 5% of the value, moving average should stay similar
    y_missing = [rand() > 0.95 ? missing : y for y in y_noisy]
    y_missing_mavg = [KolmogorovZurbenkoFilters.mavg1d(y_missing, i, 9) for i in 1:N]
    @test maximum(skipmissing(abs.(y_missing - y_missing_mavg))) < 0.2
end

@testset "kz & kza 1D" begin
    yrs = 20
    t = range(0, yrs, yrs*365)
    m = 365
    #noise
    y_error = rand(Normal(), length(t))
    y_trend = range(0,-1,length(t))
    #signal
    bkpt = 3452
    y_brk = [repeat([0.0], bkpt); repeat([0.5], length(t)-bkpt)]
    y_signal = y_trend + y_brk
    # y = seasonal + trend + break point + noise
    y = sin.(2*pi*t) + y_signal + y_error

    # kz reconstruction of signal
    y_kz = kz(y,m)
    @test maximum(abs.(y_kz - y_signal)) < maximum(abs.(y - y_signal))

    # kza reconstruction of the signal
    y_kza = kza(y,m,minimum_window_length=10, iterations=3)
    @test maximum(abs.(y_kza - y_signal)) < maximum(abs.(y - y_kz))
end