/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

Class
    Foam::mDoFSystemMotion

Description
    Multi-degree-of-freedom (mDOF) system motion for a building using 
    lamped mass system.

    Holds the mass, dampning and stifness matrix.

    The time-integrator for the motion is run-time selectable with options for
    symplectic (explicit), Crank-Nicolson and Newmark schemes.

SourceFiles
    mDoFSystemMotionI.H
    mDoFSystemMotion.C
    mDoFSystemMotionIO.C

\*---------------------------------------------------------------------------*/

#ifndef mDoFSystemMotion_H
#define mDoFSystemMotion_H

#include "mDoFSystemMotionState.H"
#include "pointField.H"
//#include "mDoFSystemMotionRestraint.H"
//#include "mDoFSystemMotionConstraint.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class mDoFSolver;

/*---------------------------------------------------------------------------*\
                      Class mDoFSystemMotion Declaration
\*---------------------------------------------------------------------------*/

class mDoFSystemMotion
{
    friend class mDoFSolver;

    // Private Data

        //- Motion state data object
        mDoFSystemMotionState motionState_;

        //- Motion state data object for previous time-step
        mDoFSystemMotionState motionState0_;

        //- Centre of mass of the initial state
        point initialCentreOfMass_;

        //- Location of origin of the building in global coordinate system.  
        point origin_;

        //- Masses of each degree of freedom (diagonal of the mass matrix).
        //  Also contains the moment of inertia
        scalarField mass_;

        //- Stifness of each degree of freedom (diagonal of the stiffness matrix).
        //  Also contains the torsional stiffness.
        scalarField stiff_;

        //- Damping of each degree of freedom (diagonal of the damping matrix).
        scalarField dampg_;

        //- Acceleration relaxation coefficient
        scalar aRelax_;

        //- Acceleration damping coefficient (for steady-state simulations)
        scalar aDamp_;

        //- Switch to turn reporting of motion data on and off
        Switch report_;


        //- Motion solver
        autoPtr<mDoFSolver> solver_;


    // Private Member Functions


        //- Update and relax accelerations from the force and torque
        void updateAcceleration(const vector& fGlobal, const vector& tauGlobal);


        // Access functions retained as private because of the risk of
        // confusion over what is a body local frame vector and what is global

        // Access

            //- Return the initial centre of rotation
            inline const point& origin() const;

        // Edit

            //- Return the centre of rotation
            inline point& origin();


public:

    //- Runtime type information
    TypeName("mDoFSystemMotion");


    // Constructors

        //- Construct null
        mDoFSystemMotion();

        //- Construct from constant and state dictionaries
        mDoFSystemMotion
        (
            const dictionary& dict,
            const dictionary& stateDict
        );

        //- Copy constructor
        mDoFSystemMotion(const mDoFSystemMotion&);


    //- Destructor
    virtual ~mDoFSystemMotion();


    // Member Functions

        // Access

            //- Return the mass
            inline scalarField mass() const;
            
            //- Return the stiffness
            inline scalarField stiff() const;

            //- Return the damping
            inline scalarField dampg() const;

            //- Return the motion state
            inline const mDoFSystemMotionState& state() const;

            //- Return the current centre of rotation
            inline const point& origin() const;

            //- Return the angular velocity in the global frame
            inline vector velocity() const;

            //- Return the current velocity
            inline const vector& v() const;

            //- Return the report Switch
            inline bool report() const;


        // Edit

            //- Store the motion state at the beginning of the time-step
            inline void newTime();



        // Constraints and Restraints


        // Update state

            //- Symplectic integration of velocities, orientation and position.
            //  Changes to Crank-Nicolson integration for subsequent iterations.
            void update
            (
                bool firstIter,
                const vector& fGlobal,
                const vector& tauGlobal,
                scalar deltaT,
                scalar deltaT0
            );

            //- Report the status of the motion
            void status() const;


        // Transformations

            //- Return the velocity of a position
            inline point displacement(const point& pt) const;

            //- Return the velocity of a position
            inline point acceleration(const point& pt) const;

                //- Return the velocity of a position
            inline point velocity(const point& pt) const;

        
            //- Transform the given initial state point by the current motion
            //  state
            inline point transform(const point& initialPoints) const;

            //- Transform the given initial state pointField by the current
            //  motion state
            tmp<pointField> transform(const pointField& initialPoints) const;

            //- Transform the given initial state pointField by the current
            //  motion state scaled by the given scale
            tmp<pointField> transform
            (
                const pointField& initialPoints,
                const scalarField& scale
            ) const;


        //- Write
        void write(Ostream&) const;

        //- Read coefficients dictionary and update system parameters,
        //  constraints and restraints but not the current state
        bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mDoFSystemMotionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
