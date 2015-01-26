def keplerian(time, p, k, ecc, omega, t0, vsys):
	from get_rvN import get_rvn
	from numpy import zeros_like

	vel = zeros_like(time)
	get_rvn(time, p, k, ecc, omega, t0, vsys, vel)
	return vel