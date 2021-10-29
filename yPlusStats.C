/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is NOT part of OpenFOAM.

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

#include "yPlusStats.H"
#include "volFields.H"
#include "turbulenceModel.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(yPlusStats, 0);
    addToRunTimeSelectionTable(functionObject, yPlusStats, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::yPlusStats::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "y+ ()");

    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "average");
    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::yPlusStats::yPlusStats
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict)
{
    read(dict);

    writeFileHeader(file());

    volScalarField* yPlusStatsPtr
    (
        new volScalarField
        (
            IOobject
            (
                typeName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        )
    );

    mesh_.objectRegistry::store(yPlusStatsPtr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::yPlusStats::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    return true;
}


bool Foam::functionObjects::yPlusStats::execute()
{
    volScalarField& yPlusStats =
        lookupObjectRef<volScalarField>(typeName);

    if (foundObject<turbulenceModel>(turbulenceModel::propertiesName))
    {
        volScalarField::Boundary& yPlusStatsBf = yPlusStats.boundaryFieldRef();

        const turbulenceModel& model =
            lookupObject<turbulenceModel>
            (
                turbulenceModel::propertiesName
            );

        const nearWallDist nwd(mesh_);
        const volScalarField::Boundary& d = nwd.y();

        // nut needed for wall function patches
        tmp<volScalarField> tnut = model.nut();
        const volScalarField::Boundary& nutBf = tnut().boundaryField();

        // U needed for plain wall patches
        const volVectorField::Boundary& UBf = model.U().boundaryField();

        const fvPatchList& patches = mesh_.boundary();

        forAll(patches, patchi)
        {
            const fvPatch& patch = patches[patchi];

            if (isA<nutWallFunctionFvPatchScalarField>(nutBf[patchi]))
            {
                const nutWallFunctionFvPatchScalarField& nutPf =
                    dynamic_cast<const nutWallFunctionFvPatchScalarField&>
                    (
                        nutBf[patchi]
                    );

                yPlusStatsBf[patchi] = nutPf.yPlus();
            }
            else if (isA<wallFvPatch>(patch))
            {
                yPlusStatsBf[patchi] =
                    d[patchi]
                   *sqrt
                    (
                        model.nuEff(patchi)
                       *mag(UBf[patchi].snGrad())
                    )/model.nu(patchi);
            }
        }
    }
    else
    {
        WarningInFunction
            << "Unable to find turbulence model in the "
            << "database: yPlusStats will not be calculated" << endl;
        return false;
    }

    return true;
}


bool Foam::functionObjects::yPlusStats::write()
{
    const volScalarField& yPlusStats =
        obr_.lookupObject<volScalarField>(type());

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << yPlusStats.name() << endl;

    const volScalarField::Boundary& yPlusStatsBf = yPlusStats.boundaryField();
    const fvPatchList& patches = mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& patch = patches[patchi];

        if (isA<wallFvPatch>(patch))
        {
            const scalarField& yPlusStatsp = yPlusStatsBf[patchi];

            const scalar minYplus = gMin(yPlusStatsp);
            const scalar maxYplus = gMax(yPlusStatsp);
            const scalar avgYplus = gAverage(yPlusStatsp);

            if (Pstream::master())
            {
                Log << "    patch " << patch.name()
                    << " y+ : min = " << minYplus << ", max = " << maxYplus
                    << ", average = " << avgYplus << nl;

                writeCurrentTime(file());
                file()
                    << token::TAB << patch.name()
                    << token::TAB << minYplus
                    << token::TAB << maxYplus
                    << token::TAB << avgYplus
                    << endl;
            }
        }
    }

    return true;
}


// ************************************************************************* //
