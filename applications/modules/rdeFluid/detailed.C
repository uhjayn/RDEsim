/*---------------------------------------------------------------------------*\
  RDEsim rdeFluid detailed chemistry
  SPDX-License-Identifier: GPL-3.0-or-later

  Species mass fractions are advected with the same Kurganov-Tadmor face
  states as rho, U and e (reconstruct(<specie>) schemes). Chemistry is a Lie
  split: transport over deltaT first, then every cell is advanced over deltaT at
  fixed density and fixed absolute internal energy by constantVolumeChemistry.
  The sensible internal energy that rdeFluid transports then changes by
  -sum_i hf_i dY_i, so the absolute energy sum_i Y_i (es_i + hf_i) is kept.
  The split is first order in deltaT.
\*---------------------------------------------------------------------------*/

#include "rdeFluid.H"
#include "fvcDiv.H"
#include "surfaceInterpolate.H"
#include "fvmDdt.H"

void Foam::solvers::rdeFluid::initialiseDetailed()
{
    multicomponent_ = dynamic_cast<multicomponentThermo*>(&thermo_);
    if (!multicomponent_)
    {
        FatalIOErrorInFunction(properties_)
            << "reactionModel detailed needs a multicomponent thermo. "
            << "constant/physicalProperties must use mixture "
            << "multicomponentMixture (psiMulticomponentThermo)."
            << exit(FatalIOError);
    }

    chemistryProperties_.set
    (
        new IOdictionary
        (
            IOobject
            (
                "chemistryProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );
    const dictionary& dict = chemistryProperties_();

    chemistryTMin_ = dict.lookupOrDefault<scalar>("TMin", 0);
    chemistryLog_ = dict.lookupOrDefault<bool>("log", false);
    fuelSpecie_ = word(dict.lookup("fuelSpecie"));
    fuelMassFractionUnburned_ = dict.lookup<scalar>("fuelMassFractionUnburned");
    if (fuelMassFractionUnburned_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "fuelMassFractionUnburned must be positive."
            << exit(FatalIOError);
    }
    if (!multicomponent_->containsSpecie(fuelSpecie_))
    {
        FatalIOErrorInFunction(dict)
            << "fuelSpecie " << fuelSpecie_ << " is not in the mechanism."
            << exit(FatalIOError);
    }

    chemistry_.set(new constantVolumeChemistry(*multicomponent_, dict));
    chemistryDeltaT_.setSize
    (
        mesh.nCells(),
        dict.lookupOrDefault<scalar>("initialChemicalTimeStep", 1e-9)
    );

    PtrList<volScalarField>& Y = multicomponent_->Y();
    forAll(Y, i)
    {
        mesh.schemes().setFluxRequired(Y[i].name());
    }

    const scalarField sumY
    (
        [&]()
        {
            scalarField sum(mesh.nCells(), 0.0);
            forAll(Y, i)
            {
                sum += Y[i].primitiveField();
            }
            return sum;
        }()
    );
    const scalar sumYMin = gMin(sumY);
    const scalar sumYMax = gMax(sumY);
    if (mag(sumYMin - 1) > 1e-6 || mag(sumYMax - 1) > 1e-6)
    {
        FatalErrorInFunction
            << "Species mass fractions do not sum to 1 in the initial fields: "
            << "min " << sumYMin << ", max " << sumYMax
            << exit(FatalError);
    }

    Info<< "rdeFluid detailed: fuel specie " << fuelSpecie_
        << ", chemistry skipped below T = " << chemistryTMin_ << " K"
        << endl;
}


void Foam::solvers::rdeFluid::speciesEquations()
{
    if (!pos.valid() || !aphiv_pos.valid() || !rho_pos.valid())
    {
        FatalErrorInFunction
            << "shockFluid has not cached the Kurganov-Tadmor face states."
            << exit(FatalError);
    }

    PtrList<volScalarField>& Y = multicomponent_->Y();
    const label defaultSpecie = multicomponent_->defaultSpecie();

    volScalarField Yt(0*Y[0]);

    forAll(Y, i)
    {
        if (i == defaultSpecie)
        {
            continue;
        }

        volScalarField& Yi = Y[i];
        const word scheme("reconstruct(" + Yi.name() + ")");

        const surfaceScalarField YPos(fvc::interpolate(Yi, pos(), scheme));
        const surfaceScalarField YNeg(fvc::interpolate(Yi, neg(), scheme));
        const surfaceScalarField phiY
        (
            "phiY_" + Yi.name(),
            aphiv_pos()*rho_pos()*YPos + aphiv_neg()*rho_neg()*YNeg
        );

        fvScalarMatrix YEqn
        (
            fvm::ddt(rho_, Yi) + fvc::div(phiY)
         ==
            fvModels().source(rho_, Yi)
        );

        YEqn.relax();
        fvConstraints().constrain(YEqn);
        YEqn.solve();
        fvConstraints().constrain(Yi);

        Yi.max(0);
        Yt += Yi;
    }

    Y[defaultSpecie] = scalar(1) - Yt;
    Y[defaultSpecie].max(0);
    multicomponent_->normaliseY();
}


void Foam::solvers::rdeFluid::chemistryStep()
{
    const scalar deltaT = runTime.deltaTValue();
    PtrList<volScalarField>& Y = multicomponent_->Y();
    volScalarField& e = thermo_.he();
    const scalarField& T = thermo.T().primitiveField();
    const scalarField& rho = rho_.primitiveField();
    const label nS = chemistry_->nSpecie();

    scalarField y(nS + 1);
    label nIntegrated = 0;
    scalar maxTDifference = 0;

    // Temperature after chemistry from the ODE, to compare with the
    // temperature the thermo recovers from the energy afterwards.
    scalarField TOde(T.size(), -1.0);

    forAll(T, celli)
    {
        if (T[celli] < chemistryTMin_)
        {
            continue;
        }

        for (label i=0; i<nS; i++)
        {
            y[i] = Y[i][celli];
        }
        y[nS] = T[celli];

        chemistry_->solve(deltaT, rho[celli], y, chemistryDeltaT_[celli]);

        scalar sum = 0;
        for (label i=0; i<nS; i++)
        {
            y[i] = max(y[i], 0.0);
            sum += y[i];
        }

        scalar dHf = 0;
        for (label i=0; i<nS; i++)
        {
            const scalar Ynew = y[i]/sum;
            dHf += chemistry_->specieThermo(i).hf()*(Ynew - Y[i][celli]);
            Y[i][celli] = Ynew;
        }

        e[celli] -= dHf;
        TOde[celli] = y[nS];
        nIntegrated++;
    }

    forAll(Y, i)
    {
        Y[i].correctBoundaryConditions();
    }

    thermo_.correct();
    pressureCorrector();

    if (chemistryLog_)
    {
        const scalarField& Tafter = thermo.T().primitiveField();
        forAll(Tafter, celli)
        {
            if (TOde[celli] > 0)
            {
                maxTDifference =
                    max(maxTDifference, mag(TOde[celli] - Tafter[celli]));
            }
        }
        Info<< "rdeFluid chemistry: " << returnReduce(nIntegrated, sumOp<label>())
            << " cells integrated, max |T(ODE) - T(thermo)| = "
            << returnReduce(maxTDifference, maxOp<scalar>()) << " K" << endl;
    }

    const volScalarField& Yfuel = multicomponent_->Y(fuelSpecie_);
    lambda_.primitiveFieldRef() = min
    (
        max(1 - Yfuel.primitiveField()/fuelMassFractionUnburned_, 0.0),
        1.0
    );
    lambda_.correctBoundaryConditions();
}


// ************************************************************************* //
