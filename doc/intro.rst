Introduction and Overview
=========================

Abstract
--------

This document starts with a brief justification for the choice of
**spatial vector algebra** (SVA) for flight simulation and
GNC (guidance–navigation and control) applications. It then presents a
complete, clean formulation of rigid-body flight dynamics using this
formalism.

All equations are written in a global **inertial** frame :math:`W`.
Therefore, Earth’s rotation influences only the atmospheric velocity,
not the equations of motion themselves. This preserves the standard
Featherstone free-body form.

The body frame is denoted :math:`B`, with its origin at the
center-of-mass (CoM) of the flying body for convenience.

The state vector is

.. math::

   x(t)
   = \left( \mathbf{p}_W(t),\; q_{WB}(t),\;
     \boldsymbol{\omega}_B(t),\;
     \mathbf{v}_B(t),\; m(t) \right),

where

- :math:`\mathbf{p}_W(t)` is the position of the CoM in the inertial frame,
- :math:`q_{WB}(t)` is the orientation quaternion of the body frame with
  respect to the inertial frame,
- :math:`\boldsymbol{\omega}_B(t)` and :math:`\mathbf{v}_B(t)` are angular
  and linear velocity vectors expressed in the body frame,
- :math:`m(t)` is the vehicle mass.


Nomenclature
------------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Symbol
     - Description
   * - :math:`W`
     - Inertial world frame (often chosen as ECI)
   * - :math:`B`
     - Body-fixed frame at the vehicle CoM
   * - :math:`I`
     - Earth–Centered Inertial (ECI) frame
   * - :math:`E`
     - Earth–fixed rotating frame (ECEF)
   * - :math:`t`
     - Time
   * - :math:`x(t)`
     - State vector :math:`(\mathbf{p}_W,q_{WB},\boldsymbol{\omega}_B,\mathbf{v}_B,m)`
   * - :math:`u(t)`
     - Control input (thrust magnitude, control torques)
   * - :math:`m(t)`
     - Vehicle mass
   * - :math:`\dot{m}_{flow}(t)`
     - Propellant mass flow rate
   * - :math:`T(t)`
     - Thrust magnitude
   * - :math:`\boldsymbol{\tau}_{ctrl,B}(t,x)`
     - Control torque in body frame
   * - :math:`\mathbf{p}_W`
     - Position of vehicle CoM in inertial frame :math:`W`
   * - :math:`\mathbf{r}_W`
     - Generic position vector expressed in :math:`W`
   * - :math:`q_{WB}`
     - Unit quaternion from body frame :math:`B` to inertial frame :math:`W`
   * - :math:`R_{WB}`
     - Rotation matrix corresponding to :math:`q_{WB}`
   * - :math:`R_{BW}`
     - Rotation matrix from :math:`W` to :math:`B`
       (:math:`R_{BW}=R_{WB}^\top`)
   * - :math:`\boldsymbol{\omega}_B`
     - Angular velocity expressed in :math:`B`
   * - :math:`\boldsymbol{\omega}_W`
     - Angular velocity expressed in :math:`W`
   * - :math:`\boldsymbol{\omega}_\oplus`
     - Earth rotation vector
   * - :math:`\mathbf{v}_B`
     - Linear velocity of CoM expressed in :math:`B`
   * - :math:`\mathbf{v}_W`
     - Linear velocity of CoM expressed in :math:`W`
   * - :math:`\mathbf{v}_{air,W}`
     - Atmospheric (wind-free) velocity in :math:`W`
   * - :math:`\mathbf{v}_{wind,W}`
     - Wind velocity in :math:`W`
   * - :math:`\mathbf{v}_{rel,W}`
     - Relative airspeed in :math:`W`
   * - :math:`\mathbf{v}`
     - Spatial motion vector (twist)
       :math:`[\boldsymbol{\omega}, \mathbf{v}]^\top`
   * - :math:`\mathbf{a}`
     - Spatial acceleration
   * - :math:`\mathbf{f}`
     - Spatial force vector (wrench)
       :math:`[\mathbf{n}, \mathbf{f}]^\top`
   * - :math:`\mathbf{v}_B`
     - Vehicle twist expressed in :math:`B`
   * - :math:`\mathbf{a}_B`
     - Vehicle spatial acceleration expressed in :math:`B`
   * - :math:`\mathbf{f}_B`
     - Net spatial wrench acting on the body in :math:`B`
   * - :math:`\mathbf{f}_{g,B}`
     - Spatial wrench due to gravity
   * - :math:`\mathbf{f}_{aero,B}`
     - Spatial wrench due to aerodynamic forces
   * - :math:`\mathbf{f}_{T,B}`
     - Spatial wrench due to thrust
   * - :math:`\mathbf{f}_{ctrl,B}`
     - Spatial wrench due to control torques
   * - :math:`\mathbf{g}_W(\mathbf{p}_W)`
     - Gravity vector expressed in :math:`W`
   * - :math:`\mathbf{F}_{g,W}, \mathbf{F}_{g,B}`
     - Gravity force in :math:`W` and :math:`B`
   * - :math:`\mathbf{F}_{aero,W}, \mathbf{F}_{aero,B}`
     - Aerodynamic force in :math:`W` and :math:`B`
   * - :math:`\mathbf{F}_{T,B}`
     - Thrust force in :math:`B`
   * - :math:`\hat{\mathbf{t}}_B`
     - Thrust unit direction in :math:`B`
   * - :math:`\rho`
     - Atmospheric density
   * - :math:`C_D`
     - Aerodynamic drag coefficient
   * - :math:`A`
     - Reference area
   * - :math:`\mathbf{I}_B(m)`
     - Spatial inertia matrix in :math:`B` (block form)
   * - :math:`\mathbf{J}_B(m)`
     - Rotational inertia matrix about the CoM in :math:`B`
   * - :math:`\mathbb{I}_F`
     - General spatial inertia in frame :math:`F`
   * - :math:`IC_F`
     - Rotational inertia at COM expressed in frame :math:`F`
   * - :math:`c_F`
     - COM location in frame :math:`F`
   * - :math:`I_3`
     - :math:`3\times 3` identity matrix
   * - :math:`[\cdot]_\times`
     - Skew-symmetric matrix for 3-D cross products
   * - :math:`\times, \times^*`
     - Spatial cross-product operator on twists and wrenches
   * - :math:`{}_A X_B`
     - Spatial motion transform from frame :math:`B` to :math:`A`
   * - :math:`{}_A X_B^*`
     - Dual force transform (co-adjoint action)
   * - :math:`\mathrm{ad}_{\mathbf{v}}`
     - Adjoint action associated with twist :math:`\mathbf{v}`
   * - :math:`\mathfrak{se}(3)`
     - Lie algebra of :math:`SE(3)` (space of twists)
   * - :math:`SE(3)`
     - Special Euclidean group (rigid-body motions)
   * - :math:`f(x,u)`
     - Nonlinear state propagation model
   * - :math:`h(x)`
     - Measurement model
   * - :math:`\mathbf{F}`
     - State-transition Jacobian :math:`\partial f/\partial x`
   * - :math:`\mathbf{H}`
     - Measurement Jacobian :math:`\partial h/\partial x`
   * - :math:`z`
     - Measurement vector
   * - :math:`v`
     - Measurement noise
   * - :math:`g_0`
     - Standard gravity
   * - :math:`I_{sp}`
     - Propellant specific impulse
   * - :math:`\omega_\oplus`
     - Earth rotation rate magnitude


Abbreviations
-------------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Abbreviation
     - Meaning
   * - AD
     - Algorithmic Differentiation
   * - CoM
     - Center of Mass
   * - DoF
     - Degree of Freedom
   * - ECEF
     - Earth-Centered Earth-Fixed (Frame :math:`E`)
   * - ECI
     - Earth-Centered Inertial (Frame :math:`I`)
   * - EKF
     - Extended Kalman Filter
   * - FDM
     - Flight Dynamics Model
   * - FMU
     - Functional Mock-up Unit
   * - GNC
     - Guidance, Navigation, and Control
   * - GPS
     - Global Positioning System
   * - IMU
     - Inertial Measurement Unit
   * - ODE
     - Ordinary Differential Equation
   * - SVA
     - Spatial Vector Algebra

