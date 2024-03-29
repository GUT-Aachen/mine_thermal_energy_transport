import pandas as pd
import pyomo.core as pyomo
from math import *


def annuity_lifetime_fact(n=30, wacc=0.03, u=None, cost_decrease=0):
    r"""
    ANNUITY WITH LIFETIME AND REPEATED INVESTMENT IMPLEMENTED ORIGINALLY IN OEMOF TOOLS ECONOMICS SUBMODULE:
    https://oemof-tools.readthedocs.io/en/stable/reference/oemof_tools.html

    Function to calculate annuity lifetime factor before input in the optimisation model, by default 0 cost_decrease is considered (investment cost remain the same)

    Parameters
    ----------
    n : int
        Horizon of the analysis, or number of years the annuity wants to be
        obtained for (n>=1)
    wacc : float
        Weighted average cost of capital (0<wacc<1)
    u : int
        Lifetime of the investigated investment. Might be smaller than the
        analysis horizon, 'n', meaning it will have to be replaced.
        Takes value 'n' if not specified otherwise (u>=1)
    cost_decrease : float
        Annual rate of cost decrease (due to, e.g., price experience curve).
        This only influences the result for investments corresponding to
        replacements, whenever u<n.
        Takes value 0, if not specified otherwise (0<cost_decrease<1)
    Returns
    -------
    float
        annuity
    """
    if u is None:
        u = n

    if ((n < 1) or (wacc < 0 or wacc > 1) or (u < 1) or
            (cost_decrease < 0 or cost_decrease > 1)):
        raise ValueError("Input arguments for 'annuity' out of bounds!")
#### annuity factor with lifetime consideration
    return ( 
        (wacc*(1+wacc)**n) / ((1 + wacc)**n - 1) *
        ((1 - ((1-cost_decrease)/(1+wacc))**n) /
         (1 - ((1-cost_decrease)/(1+wacc))**u)))

def annuity(capex, n=30, wacc=0.03, u=None, cost_decrease=0):
    r"""
    ANNUITY WITH LIFETIME AND REPEATED INVESTMENT IMPLEMENTED ORIGINALLY IN OEMOF TOOLS ECONOMICS SUBMODULE:
    https://oemof-tools.readthedocs.io/en/stable/reference/oemof_tools.html

    Calculates the annuity of an initial investment 'capex', considering
    the cost of capital 'wacc' during a project horizon 'n'

    In case of a single initial investment, the employed formula reads:

    .. math::
        \text{annuity} = \text{capex} \cdot
            \frac{(\text{wacc} \cdot (1+\text{wacc})^n)}
            {((1 + \text{wacc})^n - 1)}

    In case of repeated investments (due to replacements) at fixed intervals
    'u', the formula yields:

    .. math::
        \text{annuity} = \text{capex} \cdot
                  \frac{(\text{wacc} \cdot (1+\text{wacc})^n)}
                  {((1 + \text{wacc})^n - 1)} \cdot \left(
                  \frac{1 - \left( \frac{(1-\text{cost\_decrease})}
                  {(1+\text{wacc})} \right)^n}
                  {1 - \left(\frac{(1-\text{cost\_decrease})}{(1+\text{wacc})}
                  \right)^u} \right)

    Parameters
    ----------
    capex : float
        Capital expenditure for first investment. Net Present Value (NPV) or
        Net Present Cost (NPC) of investment
    n : int
        Horizon of the analysis, or number of years the annuity wants to be
        obtained for (n>=1)
    wacc : float
        Weighted average cost of capital (0<wacc<1)
    u : int
        Lifetime of the investigated investment. Might be smaller than the
        analysis horizon, 'n', meaning it will have to be replaced.
        Takes value 'n' if not specified otherwise (u>=1)
    cost_decrease : float
        Annual rate of cost decrease (due to, e.g., price experience curve).
        This only influences the result for investments corresponding to
        replacements, whenever u<n.
        Takes value 0, if not specified otherwise (0<cost_decrease<1)
    Returns
    -------
    float
        annuity
    """
    if u is None:
        u = n

    if ((n < 1) or (wacc < 0 or wacc > 1) or (u < 1) or
            (cost_decrease < 0 or cost_decrease > 1)):
        raise ValueError("Input arguments for 'annuity' out of bounds!")

    return (
        capex * (wacc*(1+wacc)**n) / ((1 + wacc)**n - 1) *
        ((1 - ((1-cost_decrease)/(1+wacc))**n) /
         (1 - ((1-cost_decrease)/(1+wacc))**u)))


from scipy import optimize
import numpy as np

def slope_intercept(x1,y1,x2,y2):
    a = (y2 - y1) / (x2 - x1)
    b = y1 - a * x1     
    return a,b


def segments_fit(X, Y, count):
    xmin = X.min()
    xmax = X.max()

    seg = np.full(count - 1, (xmax - xmin) / count)

    px_init = np.r_[np.r_[xmin, seg].cumsum(), xmax]
    py_init = np.array([Y[np.abs(X - x) < (xmax - xmin) * 0.01].mean() for x in px_init])

    def func(p):
        seg = p[:count - 1]
        py = p[count - 1:]
        px = np.r_[np.r_[xmin, seg].cumsum(), xmax]
        return px, py

    def err(p):
        px, py = func(p)
        Y2 = np.interp(X, px, py)
        return np.mean((Y - Y2)**2)

    r = optimize.minimize(err, x0=np.r_[seg, py_init], method='Nelder-Mead')
    return func(r.x)