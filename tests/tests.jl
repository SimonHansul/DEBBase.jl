@time using Revise
@time using DEBBase

glb = DEBBase.GlobalParams()
anm = DEBBase.DEBParams()

out = DEBBase.run_model(glb, anm)

