using Pkg; Pkg.activate("tests")
using Revise
@time using DEBBase
using Glob

tests = glob("tests/*.jl") |> 
x -> [splitpath(xi)[end] for xi in x] |>
x -> filter(f -> f != "tests.jl", x)

for test in tests
    @info("Running $test")
    include(test)
end
