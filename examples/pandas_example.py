import matplotlib.pyplot as plt
import pandas as pd

import par_est

sec = pd.date_range("31-12-1993", periods=10, freq="S")
sec_onemore = pd.date_range("31-12-1993", periods=11, freq="S")

# orig = pd.Timestamp(year=1989,month=11,day=9)
orig = pd.Timestamp(year=1970, month=1, day=1)

# time_grid = np.linspace(10, 10000, 40)
var_list, model = par_est.examples.cstr_dae()
a = var_list["e0_T"]
breakpoint()
var_res = par_est.tools.generate_varlist_with_data(var_list, model, sec.second)

times = [0, 2, 3, 4, 5.5, 8, 9, 10, 11, 12]
# from_list = pd.to_datetime(times, unit="s", origin=orig)
from_list = pd.to_datetime(times, unit="s")
print(from_list)


data = var_res.index(0).value

exp = pd.DataFrame(data, index=sec)

exp = pd.DataFrame(data, index=from_list, columns=["value"])
breakpoint()
exp.plot()
plt.show()

breakpoint()
print(sec)
