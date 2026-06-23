import mopeds


def isomerization():
    VAR_LIST, MODEL, EXP_DATA = mopeds.examples.isomerization_model()
    VAR_LIST["x1"].fixed = False
    VAR_LIST["x1"].lower_bound = 0.1
    VAR_LIST["x1"].upper_bound = 3
    VAR_LIST["theta1"].fixed = False
    oed = mopeds.OptimalExperimentalDesign_NLE(MODEL, [VAR_LIST])
    oed.optimize()
    print(oed.calculate_objective_and_jacobian({"x1": 1}))



if __name__ == "__main__":
    isomerization()
