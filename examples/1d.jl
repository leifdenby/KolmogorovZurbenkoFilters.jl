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

p1 = plot(y, label="y")
plot!(p1, y_signal, label="signal", linewidth=2)
p2 = plot(y_signal, label="signal", linewidth=2)
plot!(p2, y_kz, label="KZ applied to y", linewidth=2)
plot!(p2, y_kza, label="KZA applied to y", linewidth=2)
plot(p1, p2, layout=(2, 1), size=(800, 600))
savefig("examples/1d.png")