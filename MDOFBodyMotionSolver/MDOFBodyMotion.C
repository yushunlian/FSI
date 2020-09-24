/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2017 OpenFOAM Foundation
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

#include "MDOFBodyMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "pointPatchDist.H"
#include "pointConstraints.H"
#include "uniformDimensionedFields.H"
#include "forces.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MDOFBodyMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        MDOFBodyMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MDOFBodyMotionSolver::MDOFBodyMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    motion_
    (
        coeffDict(),
        IOobject
        (
            "MDOFBodyMotionState",
            mesh.time().timeName(),
            "uniform",
            mesh
        ).typeHeaderOk<IOdictionary>(true)
      ? IOdictionary
        (
            IOobject
            (
                "MDOFBodyMotionState",
                mesh.time().timeName(),
                "uniform",
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        )
      : coeffDict()
    ),
    patches_(wordReList(coeffDict().lookup("patches"))),
    patchSet_(mesh.boundaryMesh().patchSet(patches_)),
    di_(readScalar(coeffDict().lookup("innerDistance"))),
    do_(readScalar(coeffDict().lookup("outerDistance"))),
    test_(coeffDict().lookupOrDefault<Switch>("test", false)),
    rhoInf_(1.0),
    rhoName_(coeffDict().lookupOrDefault<word>("rho", "rho")),
    scale_
    (
        IOobject
        (
            "motionScale",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pointMesh::New(mesh),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    curTimeIndex_(-1),
    firstTime_(true)
{
    if (rhoName_ == "rhoInf")
    {
        rhoInf_ = readScalar(coeffDict().lookup("rhoInf"));
    }

    // Calculate scaling factor everywhere
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        const pointMesh& pMesh = pointMesh::New(mesh);

        pointPatchDist pDist(pMesh, patchSet_, points0());

        // Scaling: 1 up to di then linear down to 0 at do away from patches
        scale_.primitiveFieldRef() =
            min
            (
                max
                (
                    (do_ - pDist.primitiveField())/(do_ - di_),
                    scalar(0)
                ),
                scalar(1)
            );

        // Convert the scale function to a cosine
        scale_.primitiveFieldRef() =
            min
            (
                max
                (
                    0.5
                  - 0.5
                   *cos(scale_.primitiveField()
                   *Foam::constant::mathematical::pi),
                    scalar(0)
                ),
                scalar(1)
            );

        pointConstraints::New(pMesh).constrain(scale_);
        scale_.write();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::MDOFBodyMotionSolver::~MDOFBodyMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::MDOFBodyMotionSolver::curPoints() const
{
    return points0() + pointDisplacement_.primitiveField();
}


void Foam::MDOFBodyMotionSolver::solve()
{
    const Time& t = mesh().time();

    if (mesh().nPoints() != points0().size())
    {
        FatalErrorInFunction
            << "The number of points in the mesh seems to have changed." << endl
            << "In constant/polyMesh there are " << points0().size()
            << " points; in the current mesh there are " << mesh().nPoints()
            << " points." << exit(FatalError);
    }

    // Store the motion state at the beginning of the time-stepbool
    bool firstIter = false;
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        motion_.newTime();
        curTimeIndex_ = this->db().time().timeIndex();
        firstIter = true;
    }

    dimensionedVector g("g", dimAcceleration, Zero);

    if (db().foundObject<uniformDimensionedVectorField>("g"))
    {
        g = db().lookupObject<uniformDimensionedVectorField>("g");
    }
    else if (coeffDict().found("g"))
    {
        coeffDict().lookup("g") >> g;
    }

    // scalar ramp = min(max((this->db().time().value() - 5)/10, 0), 1);
    scalar ramp = 1.0;

    if (test_)
    {
        motion_.update
        (
            firstIter,
            ramp*(motion_.mass()*g.value()),
            ramp*(motion_.mass()*(motion_.momentArm() ^ g.value())),
            t.deltaTValue(),
            t.deltaT0Value()
        );
    }
    else
    {
        dictionary forcesDict;

        forcesDict.add("type", functionObjects::forces::typeName);
        forcesDict.add("patches", patches_);
        forcesDict.add("rhoInf", rhoInf_);
        forcesDict.add("rho", rhoName_);
        forcesDict.add("CofR", motion_.centreOfRotation());

        functionObjects::forces f("forces", db(), forcesDict);

        f.calcForcesMoment();

        motion_.update
        (
            firstIter,
            ramp*(f.forceEff() + motion_.mass()*g.value()),
            ramp
           *(
               f.momentEff()
             + motion_.mass()*(motion_.momentArm() ^ g.value())
            ),
            t.deltaTValue(),
            t.deltaT0Value()
        );
    	
	writeDisplacement();
    }

    // Update the displacements
    pointDisplacement_.primitiveFieldRef() =
        motion_.transform(points0(), scale_) - points0();

    // Displacement has changed. Update boundary conditions
    pointConstraints::New
    (
        pointDisplacement_.mesh()
    ).constrainDisplacement(pointDisplacement_);

}


bool Foam::MDOFBodyMotionSolver::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool valid
) const
{
    IOdictionary dict
    (
        IOobject
        (
            "MDOFBodyMotionState",
            mesh().time().timeName(),
            "uniform",
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    motion_.state().write(dict);
    return dict.regIOobject::write();
}


bool Foam::MDOFBodyMotionSolver::read()
{
    if (displacementMotionSolver::read())
    {
        motion_.read(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::MDOFBodyMotionSolver::writeDisplacement()
{
    // adjust file streams
    if (Pstream::master())
    {
	const word& fieldName = "disp";

        fileName probeDir = mesh().time().constant()/fieldName;

    	mkDir(probeDir);

        // Remove ".."
        probeDir.clean();
	

	if(firstTime_)
	{
	   // Create directory if does not exist.
	   mkDir(probeDir);
	   
	   OFstream* fPtr = new OFstream(probeDir/fieldName);

	   OFstream& fout = *fPtr;

           filePtr_.insert(fieldName, fPtr);

	   fout <<"#Time series of displacement at center of mass.\n"<< endl;

           fout<< '#' << setw(IOstream::defaultPrecision() + 6)
                << "Time" << setw(IOstream::defaultPrecision() + 7) <<"Displacement" << endl;
           
	   firstTime_ = false;
	}

	else
	{
           unsigned int w = IOstream::defaultPrecision() + 7;
           OFstream& os = *filePtr_[fieldName];

           os  << setw(w) << mesh().time().value() << setw(w) <<  motion_.centreOfMass() << endl;
	}
    }
}

// ************************************************************************* //
