using Pkg; Pkg.activate("tests")

tests = [
    raw"test01_DEB_growthrepro.jl",
    raw"test03_TD.jl"
]

for test in tests
    @info("Running $test")
    include(test)
end