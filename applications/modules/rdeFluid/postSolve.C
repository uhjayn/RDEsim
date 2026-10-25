/*---------------------------------------------------------------------------*\
  RDEsim rdeFluid objective writer
  SPDX-License-Identifier: GPL-3.0-or-later

  Writes postProcessing/rdeObjectives/objectives.dat from the objectives list
  in constant/rdeProperties. v1 accepts a scalar patch areaAverage and a
  scalar volume max, min, or average. The Python optimizer reads the same file.
\*---------------------------------------------------------------------------*/

#include "rdeFluid.H"
#include "IOstreams.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "Pstream.H"
#include "Tuple2.H"

void Foam::solvers::rdeFluid::writeObjectives() const
{
    const List<dictionary> objectives(properties_.lookup("objectives"));
    DynamicList<Tuple2<word, scalar>> rows;

    forAll(objectives, i)
    {
        const dictionary& objective = objectives[i];
        const word name(word(objective.lookup("name")));
        const word fieldName(word(objective.lookup("field")));
        const word region(word(objective.lookup("region")));
        const word operation(word(objective.lookup("operation")));

        if (!mesh.foundObject<volScalarField>(fieldName))
        {
            FatalErrorInFunction
                << "Objective " << name << " asks for field " << fieldName
                << ". v1 objective fields must be scalar volFields "
                << "(p, T, rho, lambda, or e after it has been registered)."
                << exit(FatalError);
        }
        const volScalarField& field = mesh.lookupObject<volScalarField>(fieldName);

        scalar value = 0;
        if (region == "patch")
        {
            if (operation != "areaAverage")
            {
                FatalErrorInFunction
                    << "Objective " << name << " uses patch operation " << operation
                    << ". v1 implements areaAverage."
                    << exit(FatalError);
            }
            const word patchName(word(objective.lookup("patch")));
            const label patchi = mesh.boundary().findIndex(patchName);
            if (patchi < 0)
            {
                FatalErrorInFunction
                    << "Objective " << name << " names patch " << patchName
                    << ", which is not in the mesh."
                    << exit(FatalError);
            }
            const scalarField& area = mesh.boundary()[patchi].magSf();
            const scalarField& faceValues = field.boundaryField()[patchi];
            const scalar totalArea = gSum(area);
            if (totalArea <= SMALL)
            {
                FatalErrorInFunction
                    << "Objective " << name << " patch " << patchName
                    << " has zero area."
                    << exit(FatalError);
            }
            value = gSum(faceValues*area)/totalArea;
        }
        else if (region == "volume")
        {
            if (operation == "max")
            {
                value = gMax(field.primitiveField());
            }
            else if (operation == "min")
            {
                value = gMin(field.primitiveField());
            }
            else if (operation == "average")
            {
                const scalar volume = gSum(mesh.V());
                if (volume <= SMALL)
                {
                    FatalErrorInFunction
                        << "Objective " << name << " mesh volume is zero."
                        << exit(FatalError);
                }
                value = gSum(field.primitiveField()*mesh.V())/volume;
            }
            else
            {
                FatalErrorInFunction
                    << "Objective " << name << " uses volume operation " << operation
                    << ". v1 implements max, min, and average."
                    << exit(FatalError);
            }
        }
        else
        {
            FatalErrorInFunction
                << "Objective " << name << " region '" << region
                << "' is not 'patch' or 'volume'."
                << exit(FatalError);
        }

        rows.append(Tuple2<word, scalar>(name, value));
    }

    if (!Pstream::master())
    {
        return;
    }

    OStringStream body;
    body << "# rdeFluid objectives" << nl;
    forAll(rows, i)
    {
        body << rows[i].first() << " " << rows[i].second() << nl;
    }

    const fileName directory = runTime.path()/"postProcessing"/"rdeObjectives";
    mkDir(directory);
    {
        OFstream latest(directory/"objectives.dat");
        // write(const char*) is plain text. operator<< on a Foam string quotes it.
        latest.write(body.str().c_str());
    }
    const fileName stamped = directory/runTime.name();
    mkDir(stamped);
    {
        OFstream copy(stamped/"objectives.dat");
        copy.write(body.str().c_str());
    }
}


void Foam::solvers::rdeFluid::postSolve()
{
    if (multicomponent_)
    {
        chemistryStep();
    }
    writePositivityCount();
    writeObjectives();
}


// ************************************************************************* //
