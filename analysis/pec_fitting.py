from dataclasses import dataclass
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

from potential_energy_curves.store_pecs import load_pecs, PECData


@dataclass
class ParabolaFit:
    """
    Stores parameters of parabolic fit to potential 
    well of the form:

    E(R) = E0 + a(R-R0)

    as well as the fitting region and the data that was fit.

    """
    R0 : float
    E0 : float
    a : float

    R_min : float
    R_max : float

    fitted_R : np.ndarray
    fitted_pec : np.ndarray

    R : np.ndarray
    pec : np.ndarray



def select_pec(energies, labels, state, m, channel):

    m_vals = list(energies.keys())
    if m is None:
        m = m_vals[0]

    pecs = energies[m]

    channels = list(pecs.keys())
    if channel is None:
        channel = channels[0]

    pecs = pecs[channel]

    for i, (n, l) in enumerate(labels[m][channel]):
        if (n, l) == state:
            return pecs[:, i]



def fit_parabola(
        pec_data : PECData,
        state : tuple[int, int], 
        R_min : float, 
        R_max : float, 
        m : int | None = None, 
        channel: str | None = None,
        ) -> ParabolaFit:

    R = pec_data.R
    energies = pec_data.energies
    labels = pec_data.labels

    pec = select_pec(energies, labels, state, m, channel)

    mask = (R_min < R) & (R < R_max)
    R_fit = R[mask]
    pec_fit = pec[mask]

    # Initial guesses
    coefficients = np.polyfit(R_fit, pec_fit, 2)

    a_guess = coefficients[0]
    R0_guess = -coefficients[1] / (2 * coefficients[0])
    E0_guess = np.polyval(coefficients, R0_guess)

    popt, pcov = curve_fit(
        parabola,
        R_fit,
        pec_fit,
        p0 = [E0_guess, R0_guess, a_guess]
    )

    E0, R0, a = popt

    return ParabolaFit(
        E0 = E0,
        R0 = R0,
        a = a,
        R_min = R_min,
        R_max = R_max,
        fitted_R = R_fit,
        fitted_pec = pec_fit,
        R = R,
        pec = pec
    )


def parabola(R, E0, R0, a):
    return E0 + a * (R - R0)**2

    
def plot_parabolic_fit(para_fit : ParabolaFit) -> None:

    R_para = np.linspace(para_fit.R_min, para_fit.R_max, 500)
    E_para = parabola(R_para, para_fit.E0, para_fit.R0, para_fit.a)

    label = (
        f"E0 = {para_fit.E0} \n"
        f"R0 = {para_fit.R0} \n"
        f"a = {para_fit.a}"
    )

    plt.plot(para_fit.fitted_R, para_fit.fitted_pec, "o")
    plt.plot(R_para, E_para, label = label)
    plt.grid()
    plt.legend()
    plt.show()
