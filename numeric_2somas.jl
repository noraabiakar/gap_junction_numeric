#!/usr/bin/env julia

include("HHChannels_2somas.jl")

using JSON
using Unitful
using Unitful.DefaultSymbols
using Main.HHChannels

scale(quantity, unit) = uconvert(NoUnits, quantity/unit)

radius = 18.8µm/2
area = 4*pi*radius^2
A_scale = 1µm^2

stim = Stim(20ms, 22ms, 0.1nA/area)
#conv_ggap= .007*10e-9*area/A_scale*1e-8   # expect one spike in second soma
conv_ggap=5e-5
ggap = conv_ggap*S*cm^-2

ts, vs1, vs2 = run_hh(100ms, ggap, stim=stim, sample_dt=0.025ms)

trace = Dict(
    :name => "membrane voltage",
    :sim => "numeric",
    :model => "soma",
    :units => "mV",
    :data => Dict(
        :time => scale.(ts, 1ms),
        Symbol("soma1.mid") => scale.(vs1, 1mV),
        Symbol("soma2.mid") => scale.(vs2, 1mV)
    )
)

io = open("numeric_2somas_output.txt", "w");
println(io, JSON.json([trace]))
