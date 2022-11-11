using Plots
using KolmogorovZurbenkoFilters


yrs = 20
t = range(0, yrs, yrs*365)
m = 365
# noise
y_error = rand(Normal(), length(t))
y_trend = range(0,-1,length(t))
# signal
bkpt = 3452
y_brk = [repeat([0.0], bkpt); repeat([0.5], length(t)-bkpt)]
y_signal = y_trend + y_brk
# y = seasonal + trend + break point + noise
y = sin.(2*pi*t) + y_signal + y_error

# kz reconstruction of signal
y_kz = kz(y,m)

# kza reconstruction of the signal
y_kza = kza(y,m,minimum_window_length=10, iterations=3)

plot(y, label="y", size=(1000, 400))
plot!(y_signal, label="signal", linewidth=2)
plot!(y_kz, label="KZ applied to y", linewidth=2)
plot!(y_kza, label="KZA applied to y", linewidth=2)
savefig("examples/1d.png")