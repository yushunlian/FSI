/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
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
    Foam::functionObjects::mDoFSystemState

Description
    Writes the mDoF motion state.

    Example of function object specification:
    \verbatim
    mDoFSystemState
    {
        type           mDoFSystemState;
        libs           ("libmDoFSystemState.so");
        angleFormat    degrees;
    }
    \endverbatim

Usage
    \table
        Property     | Description                  | Required | Default value
        type         | type name: mDoFSystemState    | yes |
        angleFormat  | degrees or radians           | no       | radian
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles

SourceFiles
    mDoFSystemState.C

\*---------------------------------------------------------------------------*/

#ifndef mDoFSystemState_H
#define mDoFSystemState_H

#include "fvMeshFunctionObject.H"
#include "vector.H"
#include "logFiles.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    class mDoFSystemMotion;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class mDoFSystemState Declaration
\*---------------------------------------------------------------------------*/

class mDoFSystemState
:
    public fvMeshFunctionObject,
    public logFiles
{
    // Private Data

        word angleFormat_;


    // Private Member Functions

        const mDoFSystemMotion& motion() const;


protected:

    // Protected Member Functions

        //- overloaded writeFileHeader from writeFile
        virtual void writeFileHeader(const label i = 0);


public:

    //- Runtime type information
    TypeName("mDoFSystemState");


    // Constructors

        //- Construct from Time and dictionary
        mDoFSystemState
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        mDoFSystemState(const mDoFSystemState&) = delete;


    //- Destructor
    virtual ~mDoFSystemState();


    // Member Functions

        //- Return the current body velocity
        vector velocity() const;

        //- Return the current body angular velocity
        vector angularVelocity() const;

        //- Read the mDoFSystemState data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the mDoFSystemState
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const mDoFSystemState&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
