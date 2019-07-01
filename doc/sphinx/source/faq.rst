Tips & Suggestions
------------------

**SOFT gives a 'Step-size underflow in Runge-Kutta integrator' error.**

This error is usually due to one of two reasons: (i) something is wrong with the
magnetic field you are running with, or (ii) the integrator tolerance
:option:`@Equation tolerance` is too strict.

The former reason is possible only when using a numerical magnetic field and
is indicative of a problem in one or several of the magnetic field components.
You should check to make sure that the magnetic field actually looks fine before
proceeding.

The latter reason occurs if :option:`@Equation tolerance` is set to a too small
value. Very few applications require an accuracy in the computed quantities to
machine epsilon, and for both performance and stability reasons,
:option:`@Equation tolerance` should therefore be set to a value significantly
larger than machine epsilon. On most systems, machine epsilon is around
:math:`2\cdot 10^{-16}`.
