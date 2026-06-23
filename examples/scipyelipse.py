import matplotlib.pyplot as plt
import numpy as np
import scipy
import mopeds

# v = scipy.stats.t.ppf(0.84, 10000)
# vv = scipy.stats.t.ppf(0.975, 10000)
# print(vv)
# # print(v, vv/v)
# breakpoint()

mean = [1, 0]
b = 1
cov = np.array([[2, b], [b, 1]])
rng = np.random.default_rng()
pts = rng.multivariate_normal(mean, cov, size=100000)

alpha = 0.5
# cov_scaled = 2 * scipy.stats.f(alpha, 2)
multiplier = scipy.stats.chi2.ppf(alpha, 2)
ell = mopeds.tools.CovarianceEllipse(cov * multiplier, mean, alpha)
ell.plot(x=pts[:, 0], y=pts[:, 1])
plt.show()


# breakpoint()

# ci = 0.67
# res = covariance_ellipse(cov, mean, pts[:, 0] , pts[:, 1], ci)
# fig, ax = plt.subplots()
# select = res <= 1
# print(select.sum() / select.shape[0])
# # v = ax.plot(pts[:, 0][res<1], pts[:, 1][res<=1], '.', alpha=1, c="g")
# v = ax.plot(pts[:, 0][~select], pts[:, 1][~select], '.', alpha=1, c="r")
# v = ax.plot(pts[:, 0][select], pts[:, 1][select], '.', alpha=1, c="g")
# xx, yy = plot_covariance_ellipse(cov, mean, ax, ci)

# vv = (xx.min() <= pts[:, 0]) & (pts[:, 0] <= xx.max())
# # breakpoint()

plt.show()

