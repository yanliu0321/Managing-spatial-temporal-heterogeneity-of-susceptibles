using DifferentialEquations
using LinearAlgebra, DiffEqSensitivity, Optim
using Flux: flatten, params
using DiffEqFlux, Flux
using Plots
using Flux: train!
using NNlib
using LaTeXStrings
using DataFrames
using CSV
source_data = DataFrame(CSV.File("./DeepLearningEffectiveReproductionNumber/Source_Data/SIR1.csv"))
data_acc = source_data.Total     # 取出total  这个就是每日的新增

println(length(data_acc))

n = 1
m = 129
display(plot(data_acc, label="Daily Confirmed Cases", lw=2))
savefig("./DeepLearningEffectiveReproductionNumber/Source_Data/daily.png")
trainingdata = Float32.(data_acc)                #每日确诊病例作为训练数据
# Generate Data from subexpontial growth model







ann_node = FastChain(FastDense(1, 10, tanh), FastDense(10, 1))
p0 = Float64.(initial_params(ann_node))
# println(length(p0))
# println(length([p0;[0.1]]))
#p = [p0;[0.4]]
u_0 = Float32[24350000-2.56259300554, 400002.56259300554, 1.7780148778779332, 0, 5]        #初始值和求解区间
#u_0 = Float32[24213097, 536903, 2, 0, 5]        #初始值和求解区间
tspan = (0.0f0, 129.0f0)
tsteps = range(tspan[1], tspan[2], length=length(data_acc))

p_min = [-0.1981237666251065, 0.7362071286591692, 0.1842425143035729, 0.19381538942363952, 0.2145343455544924, 0.9512717643474219, -0.7602044418381637, -0.19071589726510452, -0.02728563456495847, -1.0094412333595482, 0.36052546062529955, 0.775148779189727, -0.2594233723451719, -0.22381022804451228, -0.17773369210976697, 0.8470581355338707, -0.7763734603114396, 0.26757751152519377, 2.8213097535510583, -0.860329898812195, 0.4729363536579822, 0.5345077054263552, -0.28002305074091854, -0.2934170670257126, -0.6467468710576224, 0.5115784726534905, 0.016729394853621404, 0.29274476975528124, 1.145493103249632, -0.2137948465477729, -0.13428751894751786]


function model2_nn(du, u, p, t)
    # du[1] = 0.1 * ann_node(t, p[1:length(p0)])[1] * u[1] - p[end] * u[1]
    #du[1] = 0.1 * ann_node(t, p[1:length(p0)])[1] * u[1] - p[end] * u[1]
    S_1, S_2, I, R, D = u
    N= S_2 + I + R 
    β = 1.74
    #Enzyme.API.runtimeActivity!(true)
    du[1] = (- 1) * I * abs(ann_node(t, p)[1])
    du[2] = I * abs(ann_node(t, p)[1])- β * S_2 * I / N
    du[3] = β * S_2 * I / N - 0.2 * I
    du[4] = 0.2 * I
    du[5] = β * S_2 * I / N
end
prob_nn = ODEProblem(model2_nn, u_0, tspan, p_min)
function train(θ)
    Array(concrete_solve(prob_nn, Tsit5(), [24350000-θ[length(p0)+2];θ[length(p0)+2:end];u_0[3:4]], θ, saveat = 1,
        abstol = 1e-6, reltol = 1e-6))#,sensealg = InterpolatingAdjoint(autojacvec=ReverseDiffVJP())))
end


tspan_predict = (0.0, 129)
scatter(data_acc, label = "Real accumulated cases")
prob_prediction = ODEProblem(model2_nn, u_0, tspan, p_min)
data_prediction = Array(solve(prob_prediction, Tsit5(), saveat = 1))
pred = data_prediction[5, :]
dped = pred[(n+1):m] - pred[n:m-1]  
zpred = [pred[1]; dped] 
s1 = data_prediction[1, :]
s2 = data_prediction[2, :]
i = data_prediction[3, :]
r = data_prediction[4, :]
plt = scatter(tsteps, trainingdata, label="Accumulated cases")
plot!(plt, tsteps, zpred, label="Predicted accumulated cases")
display(plot(plt))
#savefig("./DeepLearningEffectiveReproductionNumber/Saving_Data/annepi.png")
#plot!(data_prediction[5, :], label = "Fit accumulated cases")
savefig("./DeepLearningEffectiveReproductionNumber/Saving_Data/lv.png")

println(size(data_prediction))

#println(s1)
#println(s2)
#println(i)
#println(r)

plot(s1, label="S1")
savefig("./DeepLearningEffectiveReproductionNumber/Saving_Data/s1.png")
plot(s2, label="S2")

plot!(i, label="I")

plot!(r, label="R")
savefig("./DeepLearningEffectiveReproductionNumber/Saving_Data/s2.png")
println(s1)
println(s2)
println(i)
println(r)
println(zpred)