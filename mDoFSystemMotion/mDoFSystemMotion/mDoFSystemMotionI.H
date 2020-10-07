/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::tensor Foam::mDoFSystemMotion::rotationTensorX
(
    scalar phi
) const
{
    return tensor
    (
        1, 0, 0,
        0, Foam::cos(phi), -Foam::sin(phi),
        0, Foam::sin(phi), Foam::cos(phi)
    );
}


inline Foam::tensor Foam::mDoFSystemMotion::rotationTensorY
(
    scalar phi
) const
{
    return tensor
    (
        Foam::cos(phi), 0, Foam::sin(phi),
        0, 1, 0,
        -Foam::sin(phi), 0, Foam::cos(phi)
    );
}


inline Foam::tensor Foam::mDoFSystemMotion::rotationTensorZ
(
    scalar phi
) const
{
    return tensor
    (
        Foam::cos(phi), -Foam::sin(phi), 0,
        Foam::sin(phi), Foam::cos(phi), 0,
        0, 0, 1
    );
}


inline Foam::Tuple2<Foam::tensor, Foam::vector>
Foam::mDoFSystemMotion::rotate
(
    const tensor& Q0,
    const vector& pi0,
    const scalar deltaT
) const
{
    Tuple2<tensor, vector> Qpi(Q0, pi0);
    tensor& Q = Qpi.first();
    vector& pi = Qpi.second();

    tensor R = rotationTensorX(0.5*deltaT*pi.x()/momentOfInertia_.xx());
    pi = pi & R;
    Q = Q & R;

    R = rotationTensorY(0.5*deltaT*pi.y()/momentOfInertia_.yy());
    pi = pi & R;
    Q = Q & R;

    R = rotationTensorZ(deltaT*pi.z()/momentOfInertia_.zz());
    pi = pi & R;
    Q = Q & R;

    R = rotationTensorY(0.5*deltaT*pi.y()/momentOfInertia_.yy());
    pi = pi & R;
    Q = Q & R;

    R = rotationTensorX(0.5*deltaT*pi.x()/momentOfInertia_.xx());
    pi = pi & R;
    Q = Q & R;

    return Qpi;
}


inline const Foam::PtrList<Foam::mDoFSystemMotionRestraint>&
Foam::mDoFSystemMotion::restraints() const
{
    return restraints_;
}


inline const Foam::PtrList<Foam::mDoFSystemMotionConstraint>&
Foam::mDoFSystemMotion::constraints() const
{
    return constraints_;
}


inline const Foam::point&
Foam::mDoFSystemMotion::initialCentreOfRotation() const
{
    return initialCentreOfRotation_;
}


inline const Foam::tensor&
Foam::mDoFSystemMotion::initialQ() const
{
    return initialQ_;
}


inline const Foam::tensor& Foam::mDoFSystemMotion::Q() const
{
    return motionState_.Q();
}


inline const Foam::vector& Foam::mDoFSystemMotion::v() const
{
    return motionState_.v();
}


inline const Foam::vector& Foam::mDoFSystemMotion::a() const
{
    return motionState_.a();
}


inline const Foam::vector& Foam::mDoFSystemMotion::pi() const
{
    return motionState_.pi();
}


inline const Foam::vector& Foam::mDoFSystemMotion::tau() const
{
    return motionState_.tau();
}


inline Foam::point& Foam::mDoFSystemMotion::initialCentreOfRotation()
{
    return initialCentreOfRotation_;
}


inline Foam::tensor& Foam::mDoFSystemMotion::initialQ()
{
    return initialQ_;
}


inline Foam::tensor& Foam::mDoFSystemMotion::Q()
{
    return motionState_.Q();
}


inline Foam::vector& Foam::mDoFSystemMotion::v()
{
    return motionState_.v();
}


inline Foam::vector& Foam::mDoFSystemMotion::a()
{
    return motionState_.a();
}


inline Foam::vector& Foam::mDoFSystemMotion::pi()
{
    return motionState_.pi();
}


inline Foam::vector& Foam::mDoFSystemMotion::tau()
{
    return motionState_.tau();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::scalar Foam::mDoFSystemMotion::mass() const
{
    return mass_;
}


inline const Foam::diagTensor&
Foam::mDoFSystemMotion::momentOfInertia() const
{
    return momentOfInertia_;
}


inline const Foam::mDoFSystemMotionState&
Foam::mDoFSystemMotion::state() const
{
    return motionState_;
}


inline const Foam::point& Foam::mDoFSystemMotion::centreOfRotation() const
{
    return motionState_.centreOfRotation();
}


inline const Foam::point&
Foam::mDoFSystemMotion::initialCentreOfMass() const
{
    return initialCentreOfMass_;
}


inline Foam::point Foam::mDoFSystemMotion::centreOfMass() const
{
    return transform(initialCentreOfMass_);
}


inline Foam::vector Foam::mDoFSystemMotion::momentArm() const
{
    return centreOfMass() - motionState_.centreOfRotation();
}


inline const Foam::tensor&
Foam::mDoFSystemMotion::orientation() const
{
    return Q();
}


inline Foam::vector Foam::mDoFSystemMotion::omega() const
{
    return  Q() & (inv(momentOfInertia_) & pi());
}


inline bool Foam::mDoFSystemMotion::report() const
{
    return report_;
}


inline void Foam::mDoFSystemMotion::newTime()
{
    motionState0_ = motionState_;
}


inline Foam::point& Foam::mDoFSystemMotion::centreOfRotation()
{
    return motionState_.centreOfRotation();
}


inline Foam::point Foam::mDoFSystemMotion::velocity
(
    const point& pt
) const
{
    return (omega() ^ (pt - centreOfRotation())) + v();
}


inline Foam::point Foam::mDoFSystemMotion::transform
(
    const point& initialPoint
) const
{
    return
    (
        centreOfRotation()
      + (Q() & initialQ_.T() & (initialPoint - initialCentreOfRotation_))
    );
}

// ************************************************************************* //
