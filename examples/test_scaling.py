import mopeds
import numpy as np


def dyn():
    """Tolerance of idas makes a difference. If scaling is -4, 4 doens't matter
    if -2,4 the difference in output exists
    """
    res = []
    for i in [True, False]:
        with mopeds.options(variable_scaling=i):
            vl, m = mopeds.examples.cstr(True, True, True)
            for v in vl.values():
                v.lower_bound = -2
                v.upper_bound = 1

            endtime = 1e3
            # if i:
            #     endtime = 1 / 20

            time_grid = np.linspace(0, endtime, 10)
            print("\n\n\n\n\n")
            opts = {
                "expand": 1,
                "abstol": 1e-15,
                "reltol": 1e-12,
            }
            sim = mopeds.Simulator(m, time_grid, vl, integrator_settings=opts)
            res.append(sim.simulate(algebraic=True)[2])

    print(res[1].dataframe)
    print(res[0].dataframe / res[1].dataframe)
    for key in res[0].keys():
        # print(res[0][key] / res[1][key])
        print(key)


def stead():
    res = []
    for i in [True, False]:
        with mopeds.options(variable_scaling=i):
            vl, m = mopeds.examples.cstr_nle()
            sim = mopeds.SimulatorNLE(m, vl)
            res.append(sim.simulate_fast())
    for key in res[0].keys():
        print(key)
        for j, var in enumerate(vl.get_state()):
            print(res[0][key][j, :] / res[1][key][j, :])


if __name__ == "__main__":
    dyn()
    # stead()
