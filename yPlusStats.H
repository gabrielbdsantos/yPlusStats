/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

Class
    Foam::functionObjects::yPlusStats

Group
    grpFieldFunctionObjects

Description
    Computes the near wall \f$y_1^+\f$ for turbulence models.

    Operands:
    \table
      Operand        | Type           | Location
      input          | -              | -
      output file    | dat | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/\<file\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    yPlusStats1
    {
        // Mandatory entries (unmodifiable)
        type        yPlusStats;
        libs        (yPlusStatsFunctionObject);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                             | Type | Req'd | Dflt
      type       | Type name: yPlusStats                   | word |  yes  | -
      libs       | Library name: yPlusStatsFunctionObjects | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        <solver> -postProcess -func yPlusStats
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::writeFile

SourceFiles
    yPlusStats.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_yPlusStats_H
#define functionObjects_yPlusStats_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class yPlusStats Declaration
\*---------------------------------------------------------------------------*/

class yPlusStats
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("yPlusStats");


    // Constructors

        //- Construct from Time and dictionary
        yPlusStats
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        yPlusStats(const yPlusStats&) = delete;

        //- No copy assignment
        void operator=(const yPlusStats&) = delete;


    //- Destructor
    virtual ~yPlusStats() = default;


    // Member Functions

        //- Read the yPlusStats data
        virtual bool read(const dictionary&);

        //- Calculate the yPlusStats field
        virtual bool execute();

        //- Write the yPlusStats field
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
