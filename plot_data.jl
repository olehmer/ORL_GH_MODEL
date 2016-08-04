data = readdlm("output.txt")

temps = Float64[]
ps = Float64[]

for i=1:size(data)[1]
    t_string = split(split(data[i,3],"T=")[2],",")[1]
    push!(temps,parse(Float64,t_string))

    p_string = split(split(data[i,4],"P=")[2],",")[1]
    push!(ps,parse(Float64, p_string))
end

using PyPlot

plot(temps,ps)
yscale("log")
ylim(minimum(ps),maximum(ps))
ylim(ylim()[end:-1:1])
xlabel("Temperature [K]")
ylabel("Pressure [Pa]")
title("Pure Water Vapor Atmosphere At 1 [AU] From Sun")

