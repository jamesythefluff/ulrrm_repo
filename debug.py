from physics.scattering import Scattering
from parameters import Params
from physics.wavefunctions import build_wavefunctions

params = Params()

scattering = Scattering(params.element_data.scattering_data)

wavefunctions = build_wavefunctions(params)

for R in [88.0, 88.5, 89.0]:

    k = wavefunctions[0].calculate_k(R, params)

    print(
        "R =", R,
        "k =", k,
        "a_s =", scattering.a_s(k),
        "a_p3 =", scattering.a_p(k)
    )