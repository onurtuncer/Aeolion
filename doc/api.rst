.. ------------------------------------------------------------------------------
.. Project: Aetherion
.. Copyright (c) 2025-2026, Onur Tuncer, PhD, Istanbul Technical University
..
.. SPDX-License-Identifier: MIT
.. License-Filename: LICENSE
.. ------------------------------------------------------------------------------

.. _api:

Aetherion C++ API Reference
============================

This page documents the public C++ API extracted from source headers via
Doxygen and rendered by Breathe.  All templated types are compatible with
both ``double`` and ``CppAD::AD<double>`` unless stated otherwise.

.. contents:: Modules
   :depth: 1
   :local:

Spatial Algebra
---------------

Rigid-body spatial vectors, transforms, and operators following the spatial
vector algebra conventions of Featherstone (2008).

**Storage conventions** (used throughout the library):

* :cpp:class:`Aetherion::Spatial::Twist` — packed 6-vector :math:`[\omega;\,v]`:
  indices 0–2 angular velocity [rad/s], 3–5 linear velocity [m/s].
* :cpp:class:`Aetherion::Spatial::Wrench` — packed 6-vector :math:`[M;\,F]`:
  indices 0–2 moment [N·m], 3–5 force [N] (moment-first, consistent with the
  inner product :math:`P = \mathbf{f}^\top\mathbf{v}`).

Key types and free functions:

* ``Twist<S>`` / ``Wrench<S>`` / ``Momentum<S>`` / ``Inertia<S>`` — fundamental spatial quantities.
* ``skew(v)`` — 3×3 skew-symmetric (cross-product) matrix.
* ``ad(xi)`` / ``ad_star(xi)`` / ``ad_star_times(xi, y)`` — Lie-algebra adjoint operators.
* ``CrossMotion`` / ``CrossForce`` — spatial cross-product operators.
* ``MotionTransformMatrix`` / ``ForceTransformMatrix`` / ``TransformMotion`` / ``TransformForce``
  — 6×6 spatial transforms between frames.
* ``ShiftWrenchToNewPoint`` — re-express a wrench from point P to point Q.
* ``Power`` — instantaneous spatial power :math:`P = \mathbf{f}^\top\mathbf{v}` [W].

.. doxygennamespace:: Aetherion::Spatial
   :content-only:

Lie-Group ODE Solvers
---------------------

SO(3) / SE(3) exponential maps, left Jacobians, and the RKMK family of
structure-preserving integrators.

The RKMK (Runge–Kutta–Munthe-Kaas) framework advances differential equations
that evolve on Lie groups by composing exponential maps rather than naive
Euclidean steps, preserving the manifold structure of the solution at each stage.

Sub-modules:

* ``Aetherion::ODE::RKMK::Lie`` — :math:`\mathrm{SO}(3)` and :math:`\mathrm{SE}(3)`
  exponential/logarithm maps and left Jacobians.
* ``Aetherion::ODE::RKMK::Core`` — Butcher tableaux, stage packing, Newton solver,
  and C++20 concepts constraining integrator template parameters.
* ``Aetherion::ODE::RKMK::Integrators`` — concrete Radau IIA integrators on
  :math:`SE(3)` and the product manifold :math:`SE(3) \times \mathbb{R}^n`.

**Integrator concept hierarchy** — four C++20 concepts constrain every
component of the integration stack:

.. graphviz::
   :caption: Concept hierarchy for the Aetherion RKMK integration stack

   digraph Concepts {
     graph [rankdir=TB fontname="sans-serif" fontsize=11 splines=ortho]
     node  [shape=box style="filled,rounded" fontname="sans-serif" fontsize=9]
     edge  [fontname="sans-serif" fontsize=9]

     KF  [label="«concept»\nKinematicsFieldOnSE3\<KF, Scalar\>"
          fillcolor="#D1FAE5" color="#059669"]
     VF  [label="«concept»\nVectorFieldOnProductSE3\<VF, N, Scalar\>"
          fillcolor="#FEF3C7" color="#D97706"]
     RI  [label="«concept»\nRKMKIntegratorOnProductSE3\<I, N, Scalar\>"
          fillcolor="#DBEAFE" color="#1D4ED8"]
     IF  [label="«concept»\nIntegratorFor\<I, XiField, FField, N, Scalar\>"
          fillcolor="#EDE9FE" color="#7C3AED"]
     SP  [label="SixDoFStepper\<VF, EuclidDim, IntegratorPolicy\>"
          fillcolor="#F0FDF4" color="#15803D" shape=component]
     IS  [label="ISimulator\<VF, Snapshot, EuclidDim, IntegratorPolicy\>"
          fillcolor="#F0FDF4" color="#15803D" shape=component]

     RI -> IF  [label="subsumed by" style=dashed]
     KF -> IF  [label="XiField ctor" style=dashed]
     VF -> IF  [label="FField ctor" style=dashed]
     IF -> SP  [label="requires" color="#7C3AED"]
     SP -> IS  [label="wraps"]
   }

The four concepts, in dependency order:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Concept
     - Contract
   * - ``KinematicsFieldOnSE3<KF, Scalar>``
     - ``kf(t, g, xi) → Matrix<Scalar,6,1>`` — maps body-frame twist to :math:`\mathfrak{se}(3)`
       velocity.
   * - ``VectorFieldOnProductSE3<VF, N, Scalar>``
     - ``vf(t, g, x) → Matrix<Scalar,N,1>`` — Newton-Euler right-hand side on
       :math:`SE(3) \times \mathbb{R}^N`.
   * - ``RKMKIntegratorOnProductSE3<I, N, Scalar>``
     - ``I::step(t0, g0, x0, h, opt) → StepResult``; ``StepResult`` exposes
       ``{g1, x1, converged}``.
   * - ``IntegratorFor<I, XiField, FField, N, Scalar>``
     - Combines ``RKMKIntegratorOnProductSE3`` **and**
       ``std::constructible_from<I, XiField, FField>``. This is the concept that
       gates the ``IntegratorPolicy`` template parameter of ``SixDoFStepper``
       and ``ISimulator``.

**Plugging in a custom integrator** — any type that satisfies
``IntegratorFor<I, KinematicsXiField, VF, N>`` can replace the default
Radau IIA scheme without changing any other part of the stack:

.. code-block:: cpp

   // 1. Define your integrator (must be constructible from (XiField, FField)).
   template<class XiField, class FField, int N>
   class MyExplicitRKMK4 {
   public:
       using VecE      = Eigen::Matrix<double, N, 1>;
       using StepResult = /* ... */;

       MyExplicitRKMK4(XiField xi, FField f);
       StepResult step(double t0, const SE3<double>& g0,
                       const VecE& x0, double h,
                       const NewtonOptions& opt) const;
   };

   // 2. Verify the concept at compile time (optional but recommended).
   using KFd     = RigidBody::KinematicsXiField;
   using MyVF    = RigidBody::VectorField<CentralGravityPolicy>;
   using MyIntg  = MyExplicitRKMK4<KFd, MyVF, 7>;

   static_assert(ODE::RKMK::IntegratorFor<MyIntg, KFd, MyVF, 7>);

   // 3. Drop it into the stepper or simulator.
   //    EuclidDim (second parameter) defaults to 7 and can be increased for
   //    augmented-state use cases such as Kalman-filter bias estimation.
   using MyStepper   = RigidBody::SixDoFStepper<MyVF, 7, MyIntg>;
   using MySimulator = Simulation::ISimulator<MyVF, Snapshot1, 7, MyIntg>;

   // Extended Euclidean dimension (e.g. for KF augmented states):
   // using KFStepper = RigidBody::SixDoFStepper<MyVF, 10, MyExtendedIntg>;

   // The default (Radau IIA RKMK, EuclidDim = 7) is used when all extra
   // parameters are omitted:
   using DefaultStepper   = RigidBody::SixDoFStepper<MyVF>;
   using DefaultSimulator = Simulation::ISimulator<MyVF, Snapshot1>;

.. doxygennamespace:: Aetherion::ODE::RKMK::Lie
   :content-only:

.. doxygennamespace:: Aetherion::ODE::RKMK::Core
   :content-only:

.. doxygennamespace:: Aetherion::ODE::RKMK::Integrators
   :content-only:

Flight Dynamics
---------------

State-vector construction, kinematics field, and physics policy interfaces.

**VectorField composition** — the four policy slots of
:cpp:class:`Aetherion::RigidBody::VectorField`:

.. graphviz::
   :caption: VectorField<Gravity, Aero, Thrust, Mass> — policy composition

   digraph VectorField {
     graph [rankdir=TB fontname="sans-serif" fontsize=11]
     node  [shape=box style="filled,rounded" fontname="sans-serif" fontsize=10]
     edge  [fontname="sans-serif" fontsize=9]

     VF [label="VectorField\<Gravity, Aero, Thrust, Mass\>"
         fillcolor="#D1E8FF" color="#2563EB"]

     G  [label="GravityPolicy" fillcolor="#D1FAE5" color="#059669"]
     A  [label="AeroPolicy"    fillcolor="#FEF3C7" color="#D97706"]
     T  [label="PropulsionPolicy" fillcolor="#FEF3C7" color="#D97706"]
     M  [label="MassPolicy"    fillcolor="#FCE7F3" color="#DB2777"]

     VF -> G [label="gravity"]
     VF -> A [label="aero"]
     VF -> T [label="thrust"]
     VF -> M [label="mass_model"]
   }

**Policy hierarchy** — built-in implementations for each slot:

.. graphviz::
   :caption: Gravity, Aero, Propulsion, and Mass policy families

   digraph Policies {
     graph [rankdir=LR fontname="sans-serif" fontsize=11 nodesep=0.4]
     node  [shape=box style="filled,rounded" fontname="sans-serif" fontsize=9]
     edge  [arrowhead=empty color="#555"]

     // Gravity
     GBase [label="«concept»\nGravityPolicy"  fillcolor="#D1FAE5" color="#059669"]
     G0 [label="ZeroGravityPolicy"            fillcolor="#ECFDF5" color="#059669"]
     G1 [label="CentralGravityPolicy"         fillcolor="#ECFDF5" color="#059669"]
     G2 [label="J2GravityPolicy"              fillcolor="#ECFDF5" color="#059669"]
     GBase -> G0
     GBase -> G1
     GBase -> G2

     // Aero
     ABase [label="«concept»\nAeroPolicy"     fillcolor="#FEF3C7" color="#D97706"]
     A0 [label="ZeroAeroPolicy"               fillcolor="#FFFBEB" color="#D97706"]
     A1 [label="DragOnlyAeroPolicy"           fillcolor="#FFFBEB" color="#D97706"]
     A2 [label="BrickDampingAeroPolicy"       fillcolor="#FFFBEB" color="#D97706"]
     A3 [label="WindAwareDragPolicy\<Wind\>"  fillcolor="#FFFBEB" color="#D97706"]
     ABase -> A0
     ABase -> A1
     ABase -> A2
     ABase -> A3

     // Wind models (sub-hierarchy)
     WBase [label="is_wind_model\<W\>" fillcolor="#FFF7ED" color="#EA580C"]
     W0 [label="ZeroWind"             fillcolor="#FFF7ED" color="#EA580C"]
     W1 [label="ConstantECEFWind"     fillcolor="#FFF7ED" color="#EA580C"]
     W2 [label="PowerLawWindShear"    fillcolor="#FFF7ED" color="#EA580C"]
     W3 [label="GeodesicCallbackWind" fillcolor="#FFF7ED" color="#EA580C"]
     WBase -> W0; WBase -> W1; WBase -> W2; WBase -> W3
     A3 -> WBase [label="Wind" style=dashed arrowhead=open]

     // Propulsion
     TBase [label="«concept»\nPropulsionPolicy" fillcolor="#EDE9FE" color="#7C3AED"]
     T0 [label="ZeroPropulsionPolicy"           fillcolor="#F5F3FF" color="#7C3AED"]
     T1 [label="ConstantThrustPolicy"           fillcolor="#F5F3FF" color="#7C3AED"]
     TBase -> T0; TBase -> T1

     // Mass
     MBase [label="«concept»\nMassPolicy"  fillcolor="#FCE7F3" color="#DB2777"]
     M0 [label="ConstantMassPolicy"        fillcolor="#FDF2F8" color="#DB2777"]
     M1 [label="LinearBurnPolicy"          fillcolor="#FDF2F8" color="#DB2777"]
     MBase -> M0; MBase -> M1
   }

**Policy system** — :cpp:class:`Aetherion::RigidBody::VectorField` is templated on
four policy types.  Each policy must satisfy the corresponding C++20 concept:

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Concept
     - Required signature
     - Built-in policies
   * - ``GravityPolicy``
     - ``Wrench<S> operator()(const SE3<S>&, S mass) const``
     - ``ZeroGravityPolicy``, ``CentralGravityPolicy``, ``J2GravityPolicy``
   * - ``AeroPolicy``
     - ``Wrench<S> operator()(const SE3<S>&, const Matrix<S,6,1>&, S mass, S t) const``
     - ``ZeroAeroPolicy``
   * - ``PropulsionPolicy``
     - same as ``AeroPolicy``
     - ``ZeroPropulsionPolicy``, ``ConstantThrustPolicy``
   * - ``MassPolicy``
     - ``S mdot(S t, S mass) const``
     - ``ConstantMassPolicy``, ``LinearBurnPolicy``

All policies satisfy their concept for both ``S = double`` and
``S = CppAD::AD<double>`` (enforced by ``static_assert`` at definition time).

:cpp:class:`Aetherion::RigidBody::KinematicsXiField` implements the
right-trivialised kinematic equation
:math:`\dot{g} = g\,\hat{\xi}` as a stateless identity pass-through
(the body-frame twist *is* the Lie-algebra velocity in this formulation).

.. doxygennamespace:: Aetherion::FlightDynamics
   :content-only:

Rigid Body
----------

The 6-DoF rigid-body simulation layer built on top of the Spatial and
FlightDynamics modules.

**State manifold** :math:`SE(3) \times \mathbb{R}^7`:

* ``g ∈ SE(3)`` — attitude (rotation matrix) and ECI position.
* ``ν_B ∈ ℝ⁶`` — body-frame twist :math:`[\omega_B;\,v_B]` [rad/s | m/s].
* ``m ∈ ℝ`` — current vehicle mass [kg].

Key types:

* :cpp:struct:`Aetherion::RigidBody::State` — full state on :math:`SE(3) \times \mathbb{R}^7`.
* :cpp:struct:`Aetherion::RigidBody::InertialParameters` — mass, inertia tensor about CoG, CoG offset (all in body frame).
* :cpp:class:`Aetherion::RigidBody::VectorField` — Newton-Euler ODE right-hand side, templated on physics policies.
* :cpp:class:`Aetherion::RigidBody::SixDoFStepper` — high-level facade with three template parameters:
  ``VectorField``, ``EuclidDim`` (default 7), and ``IntegratorPolicy`` (default Radau IIA RKMK),
  all constrained by ``IntegratorFor``.
* :cpp:struct:`Aetherion::RigidBody::StateLayout` — flat 14-element index map for serialised state vectors.

.. note::

   ``SixDoFStepper<VF>`` is a shorthand for
   ``SixDoFStepper<VF, 7, RadauIIA_RKMK_ProductSE3<KinematicsXiField, VF, 7>>``.
   Pass an explicit ``EuclidDim`` (≥ 7) to extend the Euclidean state for augmented-state
   applications (e.g. Kalman-filter bias states occupy slots 7…EuclidDim−1;
   ``pack()`` zero-initialises them and ``unpack()`` ignores them).
   Pass an explicit ``IntegratorPolicy`` to swap in any integrator satisfying
   ``IntegratorFor<I, KinematicsXiField, VF, EuclidDim>`` — see the
   :ref:`Lie-Group ODE Solvers <api>` section for a worked example.

**State vector representations**

Two complementary representations are used throughout the library.

*Manifold representation* — :cpp:struct:`Aetherion::RigidBody::State`

The state lives on the product manifold :math:`SE(3) \times \mathbb{R}^7`:

.. math::

   x \;=\; \bigl(g,\;\nu_B,\;m\bigr)
   \;\in\; SE(3) \times \mathbb{R}^6 \times \mathbb{R},

where

- :math:`g = (R_{WB},\,\mathbf{p}_W) \in SE(3)` — pose: body-to-inertial rotation matrix
  and ECI position of the CoM.
- :math:`\nu_B = [\boldsymbol{\omega}_B;\,\mathbf{v}_B] \in \mathbb{R}^6` — body-frame
  twist: angular velocity [rad/s] stacked above linear velocity [m/s].
- :math:`m \in \mathbb{R}` — current vehicle mass [kg].

The Lie-group part :math:`g` is advanced by the RKMK integrator directly on the
:math:`SE(3)` manifold; the Euclidean part :math:`[\nu_B;\,m] \in \mathbb{R}^7` is
advanced by the standard Euclidean ODE solver.

*Serialised representation* — :cpp:struct:`Aetherion::RigidBody::StateLayout`

For file I/O, finite-difference Jacobians, or any context that requires a
contiguous :math:`\mathbb{R}^{14}` vector, the state is packed in this order:

.. list-table::
   :header-rows: 1
   :widths: 12 18 10 60

   * - Indices
     - Constant
     - Size
     - Quantity
   * - 0–2
     - ``IDX_P``
     - 3
     - ECI position :math:`\mathbf{p}_W` [m]
   * - 3–6
     - ``IDX_Q``
     - 4
     - Unit quaternion :math:`q_{WB}` stored as :math:`(w,\,x,\,y,\,z)` — scalar part first
   * - 7–9
     - ``IDX_W``
     - 3
     - Body angular velocity :math:`\boldsymbol{\omega}_B` [rad/s]
   * - 10–12
     - ``IDX_V``
     - 3
     - Body linear velocity :math:`\mathbf{v}_B` [m/s]
   * - 13
     - ``IDX_M``
     - 1
     - Vehicle mass :math:`m` [kg]

Total: :math:`N = 14` elements (``StateLayout::N``).

.. note::

   The quaternion :math:`q_{WB}` represents the rotation that maps vectors from body
   frame :math:`B` to inertial frame :math:`W`.  Storage order is scalar-first:
   :math:`(w, x, y, z)`, satisfying :math:`w^2 + x^2 + y^2 + z^2 = 1`.

   The corresponding rotation matrix is :math:`R_{WB}` stored in :math:`g`.
   Both :math:`g` (manifold) and the serialised :math:`q_{WB}` (flat vector)
   encode the same orientation; ``BuildInitialState`` populates them from the
   same configuration, and they remain consistent throughout integration.

.. doxygennamespace:: Aetherion::RigidBody
   :content-only:

Environment
-----------

Atmospheric models, gravity acceleration models, and WGS-84 geodetic constants.

* **US Standard Atmosphere 1976** — density, pressure, temperature, and speed of
  sound as functions of altitude.  CppAD-friendly.
* **Gravity models** — ``CentralGravity(r, mu)`` (inverse-square) and
  ``J2(r, mu, Re, J2)`` (first zonal harmonic).
* **WGS-84 constants** — ``kSemiMajorAxis_m``, ``kFlattening``, ``kGM_m3_s2``,
  ``kJ2``, ``kEarthAngularVelocity_rad_s`` and others in ``Aetherion::Environment::WGS84``.

.. doxygennamespace:: Aetherion::Environment
   :content-only:

Coordinate Transforms
---------------------

Geodetic ↔ ECEF ↔ ECI conversions, NED frame utilities, and attitude decomposition.

**Frame definitions** used throughout:

* **ECI** (Earth-Centred Inertial) — origin at Earth's centre, axes fixed to
  distant stars; used as the integration frame.
* **ECEF** (Earth-Centred Earth-Fixed) — same origin, rotates with Earth at
  :math:`\dot{\theta}_{ERA}`.
* **NED** (North-East-Down) — local tangent frame at a surface point;
  right-handed with x = North, y = East, z = Down.
* **Body** — x forward, y right, z down.

All functions accept a template scalar ``S`` and are CppAD-friendly.

**Transform chain** (forward — local to inertial):

.. code-block:: text

   Geodetic (φ, λ, h)
       → ECEF  [GeodeticToECEF]
       → ECI   [ECEFToECI, requires θ_ERA]
       → launch state [MakeLaunchStateECI]

**Transform chain** (inverse — inertial to local):

.. code-block:: text

   ECI position
       → ECEF  [ECIToECEF]
       → Geodetic (φ, λ, h)  [ECEFToGeodeticWGS84]
       → NED   [ECEFToNED / ECIToNED]
       → local orientation  [QuaternionToAzZenRollNED]

.. doxygennamespace:: Aetherion::Coordinate
   :content-only:

Serialization
-------------

JSON (de)serialization for all configuration and state types, built on
`nlohmann/json <https://github.com/nlohmann/json>`_.

* ``LoadConfig(filename)`` / ``LoadConfig(json)`` — load a
  :cpp:struct:`Aetherion::RigidBody::Config` from a JSON file or a pre-parsed
  ``nlohmann::json`` object.
* ``LoadConfigFromString(str)`` — parse and validate a JSON string directly
  (throws ``std::runtime_error`` on empty input, parse failure, or schema mismatch).
* ``from_json`` / ``to_json`` specialisations for all configuration structs:
  ``GeodeticPoseNED``, ``VelocityNED``, ``BodyRates``, ``InertialParameters``,
  ``AerodynamicParameters``, and ``Config``.

.. doxygennamespace:: Aetherion::Serialization
   :content-only:

Simulation Framework
--------------------

Application base class, argument parser, simulator interface, and telemetry types.

**Simulator class hierarchy** — all concrete simulators inherit from
:cpp:class:`Aetherion::Simulation::ISimulator`:

The full signature is ``ISimulator<VF, Snapshot, EuclidDim, IntegratorPolicy>`` where
``EuclidDim`` defaults to 7 and ``IntegratorPolicy`` defaults to ``RadauIIA_RKMK_ProductSE3``.
Existing two-parameter instantiations ``ISimulator<VF, Snapshot>`` are unchanged.
Increase ``EuclidDim`` to extend the Euclidean state (e.g. for Kalman-filter augmentation).
Pass a custom ``IntegratorPolicy`` to use an alternative integrator.

``TwoStageRocketSimulator<EuclidDim, IntegratorPolicy>`` uses the same parameters
directly (without inheriting ``ISimulator``) because it requires mutable access
to the state for the discontinuous staging mass drop.

.. graphviz::
   :caption: ISimulator inheritance — one concrete class per NASA scenario group

   digraph Simulators {
     graph [rankdir=TB fontname="sans-serif" fontsize=11]
     node  [shape=box style="filled,rounded" fontname="sans-serif" fontsize=9]
     edge  [arrowhead=empty color="#555"]

     IS [label="ISimulator\<VF, Snapshot, EuclidDim, IntegratorPolicy\>"
         fillcolor="#DBEAFE" color="#1D4ED8" shape=box]

     DS  [label="DraglessSphereSimulator\n(Scenarios 1, 9, 10)"
          fillcolor="#EFF6FF" color="#1D4ED8"]
     SW  [label="SphereWithAtmosphericDragSimulator\n(Scenario 6)"
          fillcolor="#EFF6FF" color="#1D4ED8"]
     TN  [label="TumblingBrickNoDampingSimulator\n(Scenario 2)"
          fillcolor="#EFF6FF" color="#1D4ED8"]
     TW  [label="TumblingBrickWithDampingSimulator\n(Scenario 3)"
          fillcolor="#EFF6FF" color="#1D4ED8"]
     SSW [label="DroppedSphereSteadyWindSimulator\n(Scenario 7)"
          fillcolor="#EFF6FF" color="#1D4ED8"]
     S2W [label="DroppedSphere2DWindShearSimulator\n(Scenario 8)"
          fillcolor="#EFF6FF" color="#1D4ED8"]

     IS -> DS
     IS -> SW
     IS -> TN
     IS -> TW
     IS -> SSW
     IS -> S2W
   }

**Application class hierarchy** — concrete applications follow the
**Template Method** pattern defined in
:cpp:class:`Aetherion::Simulation::Application`:

.. graphviz::
   :caption: Application inheritance — one concrete subclass per scenario executable

   digraph Applications {
     graph [rankdir=TB fontname="sans-serif" fontsize=11]
     node  [shape=box style="filled,rounded" fontname="sans-serif" fontsize=9]
     edge  [arrowhead=empty color="#555"]

     App [label="Application\n(defines run() loop + virtual hooks)"
          fillcolor="#D1FAE5" color="#065F46" shape=box]

     DA  [label="DraglessSphereApplication\n(Scenario 1)"
          fillcolor="#ECFDF5" color="#065F46"]
     SWA [label="SphereWithAtmosphericDragApplication\n(Scenario 6)"
          fillcolor="#ECFDF5" color="#065F46"]
     TNA [label="TumblingBrickNoDampingApplication\n(Scenario 2)"
          fillcolor="#ECFDF5" color="#065F46"]
     TWA [label="TumblingBrickWithDampingApplication\n(Scenario 3)"
          fillcolor="#ECFDF5" color="#065F46"]
     ECA [label="EastwardCannonballApplication\n(Scenario 9)"
          fillcolor="#ECFDF5" color="#065F46"]
     NCA [label="NorthwardCannonballApplication\n(Scenario 10)"
          fillcolor="#ECFDF5" color="#065F46"]
     SWW [label="DroppedSphereSteadyWindApplication\n(Scenario 7)"
          fillcolor="#ECFDF5" color="#065F46"]
     S2A [label="DroppedSphere2DWindShearApplication\n(Scenario 8)"
          fillcolor="#ECFDF5" color="#065F46"]

     App -> DA;  App -> SWA; App -> TNA; App -> TWA
     App -> ECA; App -> NCA; App -> SWW; App -> S2A
   }

**Template Method pattern** — :cpp:class:`Aetherion::Simulation::Application` defines
the simulation loop in ``run()`` and exposes the following virtual hooks for
subclasses to override:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Hook
     - Responsibility
   * - ``prepareSimulation()``
     - Allocate and initialise all simulation objects before the loop starts.
   * - ``writeInitialSnapshot(ofstream&)``
     - Write the :math:`t=0` state to the output CSV.
   * - ``stepAndRecord(ofstream&, double t, bool doWrite)``
     - Advance by one time step; write to CSV only when ``doWrite`` is true.
   * - ``logFinalSummary()``
     - Print end-of-simulation diagnostics.

:cpp:class:`Aetherion::Simulation::ArgumentParser` parses ``--timeStep``,
``--startTime``, ``--endTime``, ``--inputFileName``, ``--outputFileName``, and
``--writeInterval`` from ``argv``, throwing ``std::invalid_argument`` on
unknown flags, missing values, or constraint violations.

.. doxygennamespace:: Aetherion::Simulation
   :content-only:
