"""Extended functionality for GEKKO. Hopefully will be obsolete if these
extensions get implemented in the core GEKKO system.

Usage
-----
Register all extensions (floor, switch, max) with

        m = register_extensions(GEKKO())
"""
import math
from functools import partial


def register_extensions(m):
    """Registers all extensions defined in this file as GEKKO functions.

    Parameters
    ----------
    m : GEKKO
        The gekko model.

    Returns
    -------
    m : GEKKO
        The gekko model with the extensions added as functions.
    """
    m.floor = partial(_floor, m=m)
    m.ceil = partial(_ceil, m=m)
    m.switch = partial(_switch, m=m)
    m.max = partial(_max, m=m)
    m.min = partial(_min, m=m)
    m.clamp = partial(_clamp, m=m)

    return m


def _floor(x, n=10, m=None):
    """Approximates the floor function using gekko.

    Parameters
    ----------
    x : gekko variable or intermediate
        The value of which to take the floor.
    n : int
        The precision of the floor operation. The larger, the more precise.
    m : GEKKO
        A gekko model.

    Returns
    -------
    fx : gekko intermediate
        A gekko intermediate with an approximation of floor(x).
    """
    if m is None:
        raise Exception('Missing GEKKO model definition.')

    intermediates = [x - 0.5]
    pi = math.pi
    for k in range(1, n):
        intermediates.append((1 / pi) * m.sin(2 * pi * k * x) / k)

    return m.Intermediate(sum(intermediates))


def _ceil(x, n=10, m=None):
    """Approximates the ceil function using gekko.

    Parameters
    ----------
    x : gekko variable or intermediate
        The value of which to take the floor.
    n : int
        The precision of the floor operation. The larger, the more precise.
    m : GEKKO
        A gekko model.

    Returns
    -------
    fx : gekko intermediate
        A gekko intermediate with an approximation of floor(x).
    """
    return m.Intermediate(-1 * _floor(-1 * x, n, m))


def _switch(left, right, on, loc, k=100, m=None):
    """Approximates the switch function using gekko.

    The switch function returns left if `on <= loc`, and right otherwise.

    Parameters
    ----------
    left : gekko variable or intermediate
    right : gekko variable or intermediate
    on : gekko variable or intermediate
    loc : gekko variable or intermediate
    k : int >= 1, default = 100
        The quality of the switch approximation. Larger is better.
    m : GEKKO

    Returns
    -------
    fx : gekko intermediate
        An approximation of the switch function.
    """
    if m is None:
        raise Exception('Missing GEKKO model definition.')

    sig = 1 / (1 + m.exp(-k * (on - loc)))
    return m.Intermediate((1 - sig) * left + sig * right)


def _clamp(f, lower, upper, k=100, m=None):
    """Clamps the value of f such that the result is lower if f < lower, and
    upper if f > upper, or f otherwise.

    Parameters
    ----------
    f : gekko variable or intermediate
    lower : gekko variable or intermediate or number
    upper : gekko variable or intermediate or number
    k : int >= 1, default = 100
        The quality of the switch approximation. Larger is better.
    m : GEKKO

    returns
    -------
    fx : gekko intermediate
        An approximation of the clamp function.
    """
    if m is None:
        raise Exception('Missing GEKKO model definition.')

    fup = m.switch(f, upper, f, upper)
    return m.switch(lower, fup, fup, lower)


def _max(left, right, k=100, m=None):
    """Approximates the max function using gekko.

    The `max(left, right)` returns `right` if `left <= right`, and `left`
    otherwise.

    Parameters
    ----------
    left : gekko variable or intermediate
    right : gekko variable or intermediate
    k : int > 100
        The quality of the switch approximation with which max is implemented.
        Larger is better.
    m : GEKKO

    Returns
    -------
    fx : gekko intermediate
        An approximation of the max function.
    """
    if m is None:
        raise Exception('Missing GEKKO model definition.')

    return m.Intermediate(_switch(right, left, left, right, k, m))


def _min(left, right, k=100, m=None):
    """Approximates the min function using gekko.

    The `max(left, right)` returns `left` if `left <= right`, and `right`
    otherwise.

    Parameters
    ----------
    left : gekko variable or intermediate
    right : gekko variable or intermediate
    k : int > 100
        The quality of the switch approximation with which nub is implemented.
        Larger is better.
    m : GEKKO

    Returns
    -------
    fx : gekko intermediate
        An approximation of the min function.
    """
    if m is None:
        raise Exception('Missing GEKKO model definition.')

    return m.Intermediate(_switch(left, right, left, right, k, m))
