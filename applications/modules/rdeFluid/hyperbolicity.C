/*---------------------------------------------------------------------------*\
  RDEsim rdeFluid hyperbolicity check
  SPDX-License-Identifier: GPL-3.0-or-later

  The reactive Euler equations with a perfect gas are hyperbolic when gamma,
  density, pressure, and c^2 stay positive and the source is local. Floors and
  the allowed model lists come from constant/rdeProperties. A failure is written
  to postProcessing/hyperbolicity/report and then aborts the run.
\*---------------------------------------------------------------------------*/

#include "rdeFluid.H"
#include "IOdictionary.H"
#include "IOstreams.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "Pstream.H"

namespace
{

bool listed(const Foam::wordList& words, const Foam::word& key)
{
    forAll(words, i)
    {
        if (words[i] == key)
        {
            return true;
        }
    }
    return false;
}

Foam::string join(const Foam::wordList& words)
{
    Foam::string text;
    forAll(words, i)
    {
        if (i)
        {
            text += ", ";
        }
        text += words[i];
    }
    return text;
}

} // namespace


void Foam::solvers::rdeFluid::checkHyperbolicity() const
{
    const dictionary& hyper = properties_.subDict("hyperbolicity");
    const bool strict = hyper.lookup<bool>("requireStrictlyHyperbolic");
    const scalar gammaFloor = hyper.lookup<scalar>("gammaMin");
    const scalar rhoFloor = hyper.lookup<scalar>("rhoMin");
    const scalar pFloor = hyper.lookup<scalar>("pMin");
    const scalar tFloor = hyper.lookup<scalar>("TMin");
    const scalar c2Floor = hyper.lookup<scalar>("c2Min");
    const scalar maxCoStability = hyper.lookup<scalar>("maxCoStability");
    const bool sourceOnGradients =
        hyper.lookup<bool>("sourceDependsOnGradients");
    const word coupling(word(hyper.lookup("boundaryCoupling")));
    const wordList fluxAllowed(hyper.lookup("hyperbolicFluxSchemes"));
    const wordList eosAllowed(hyper.lookup("hyperbolicEquationsOfState"));
    const wordList sourceAllowed(hyper.lookup("localSourceIntegrations"));

    wordList fluxImplemented(2);
    fluxImplemented[0] = "Kurganov";
    fluxImplemented[1] = "Tadmor";
    wordList sourceImplemented(2);
    sourceImplemented[0] = "pointImplicit";
    sourceImplemented[1] = "explicit";
    wordList eosImplemented(1);
    eosImplemented[0] = "perfectGas";

    const volScalarField gammaField(thermo.Cp()/thermo.Cv());
    const scalar gammaMinField = gMin(gammaField.primitiveField());
    const scalar rhoMinField = gMin(rho_.primitiveField());
    const scalar pMinField = gMin(p_.primitiveField());
    const scalar tMinField = gMin(thermo.T().primitiveField());

    scalar c2MinField = 0;
    const bool psiPositive = gMin(thermo.psi().primitiveField()) > 0;
    if (psiPositive)
    {
        const volScalarField c2Field(gammaField/thermo.psi());
        c2MinField = gMin(c2Field.primitiveField());
    }

    DynamicList<string> failures;
    DynamicList<string> warnings;

    if (gammaMinField <= gammaFloor)
    {
        OStringStream message;
        message
            << "gamma " << gammaMinField
            << " is not above hyperbolicity.gammaMin " << gammaFloor;
        failures.append(message.str());
    }
    if (rhoMinField <= rhoFloor)
    {
        failures.append("Density is not above hyperbolicity.rhoMin");
    }
    if (pMinField <= pFloor)
    {
        failures.append("Pressure is not above hyperbolicity.pMin");
    }
    if (tMinField <= tFloor)
    {
        failures.append("Temperature is not above hyperbolicity.TMin");
    }
    if (!psiPositive || c2MinField <= c2Floor)
    {
        failures.append
        (
            "c^2 is not above hyperbolicity.c2Min, so the flux Jacobian is not real"
        );
    }

    if (!listed(eosAllowed, equationOfState_))
    {
        failures.append
        (
            "Equation of state '" + equationOfState_
          + "' is not in hyperbolicity.hyperbolicEquationsOfState"
        );
    }
    if (!listed(eosImplemented, equationOfState_))
    {
        failures.append
        (
            "Equation of state '" + equationOfState_
          + "' is not implemented. This version writes " + join(eosImplemented)
          + "."
        );
    }
    if (thermo.thermoName().find(equationOfState_) == string::npos)
    {
        failures.append
        (
            "thermo.thermoName() does not contain equationOfState '"
          + equationOfState_ + "'. The runtime thermo is " + thermo.thermoName()
        );
    }

    if (!listed(fluxAllowed, fluxScheme))
    {
        failures.append
        (
            "Flux scheme '" + fluxScheme
          + "' is not in hyperbolicity.hyperbolicFluxSchemes"
        );
    }
    if (!listed(fluxImplemented, fluxScheme))
    {
        failures.append
        (
            "Flux scheme '" + fluxScheme
          + "' is not implemented by rdeFluid. Implemented schemes: "
          + join(fluxImplemented) + "."
        );
    }

    if (!listed(sourceAllowed, sourceIntegration_))
    {
        failures.append
        (
            "Source integration '" + sourceIntegration_
          + "' is not in hyperbolicity.localSourceIntegrations"
        );
    }
    if (!listed(sourceImplemented, sourceIntegration_))
    {
        failures.append
        (
            "Source integration '" + sourceIntegration_
          + "' is not implemented. Implemented: " + join(sourceImplemented) + "."
        );
    }

    if (sourceOnGradients)
    {
        failures.append
        (
            "hyperbolicity.sourceDependsOnGradients is yes. A source that "
            "depends on spatial derivatives is not a hyperbolic balance law."
        );
    }
    if (coupling != "local")
    {
        failures.append
        (
            "hyperbolicity.boundaryCoupling is '" + coupling
          + "'. Only 'local' keeps the boundary treatment from adding an "
            "elliptic coupling."
        );
    }

    if (flowModel_ != "euler" && flowModel_ != "navierStokes")
    {
        failures.append("Unknown flowModel '" + flowModel_ + "'");
    }
    else if (flowModel_ == "navierStokes" || !inviscid)
    {
        if (strict)
        {
            failures.append
            (
                "Viscous transport is active and "
                "hyperbolicity.requireStrictlyHyperbolic is yes. "
                "Viscous terms make the system incompletely parabolic. "
                "Set flowModel euler with zero viscosity, or set "
                "requireStrictlyHyperbolic no."
            );
        }
        else
        {
            warnings.append
            (
                "Navier-Stokes transport is incompletely parabolic. "
                "The inviscid flux is still the reactive Euler flux."
            );
        }
    }

    if (ddtScheme_ == "steadyState" || steady)
    {
        if (strict)
        {
            failures.append
            (
                "ddtScheme steadyState is not a hyperbolic time evolution"
            );
        }
        else
        {
            warnings.append
            (
                "ddtScheme steadyState does not advance a time-hyperbolic system"
            );
        }
    }

    if (runTime.controlDict().found("maxCo"))
    {
        const scalar maxCo = runTime.controlDict().lookup<scalar>("maxCo");
        if (maxCo > maxCoStability)
        {
            OStringStream message;
            message
                << "controlDict maxCo " << maxCo
                << " is above hyperbolicity.maxCoStability " << maxCoStability
                << ". That is a stability limit for the explicit flux, "
                << "not a loss of hyperbolicity.";
            warnings.append(message.str());
        }
    }
    else
    {
        warnings.append
        (
            "controlDict has no maxCo, so the CFL stability comparison was not made."
        );
    }

    typeIOobject<IOdictionary> acousticHeader
    (
        "acousticReflection",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );
    if (acousticHeader.headerOk())
    {
        IOdictionary acoustic(acousticHeader);
        if (acoustic.found("timeDomainModel"))
        {
            const word model(word(acoustic.lookup("timeDomainModel")));
            if (model != "slip" && model != "resistive" && model != "helmholtz")
            {
                failures.append
                (
                    "Acoustic time-domain model '" + model + "' is not local"
                );
            }
        }
    }

    const bool hyperbolic = failures.empty();

    OStringStream report;
    report << (hyperbolic ? "Hyperbolicity: hyperbolic" : "Hyperbolicity: NOT hyperbolic")
        << nl
        << "  gamma                 " << gammaMinField << nl
        << "  density               " << rhoMinField << " kg/m^3" << nl
        << "  pressure              " << pMinField << " Pa" << nl
        << "  temperature           " << tMinField << " K" << nl
        << "  c^2                   " << c2MinField << " m^2/s^2" << nl
        << "  flux scheme           " << fluxScheme << nl
        << "  flow model            " << flowModel_ << nl
        << "  equation of state     " << equationOfState_ << nl
        << "  thermo                " << thermo.thermoName() << nl;
    if (failures.size())
    {
        report << "Failures:" << nl;
        forAll(failures, i)
        {
            report << "  - " << failures[i] << nl;
        }
    }
    if (warnings.size())
    {
        report << "Warnings:" << nl;
        forAll(warnings, i)
        {
            report << "  - " << warnings[i] << nl;
        }
    }
    if (!failures.size() && !warnings.size())
    {
        report
            << "The selected equation of state, source, and scheme satisfy the checks."
            << nl;
    }

    if (Pstream::master())
    {
        const fileName directory =
            runTime.path()/"postProcessing"/"hyperbolicity";
        mkDir(directory);
        OFstream os(directory/"report");
        os << report.str();
    }

    Info<< report.str() << endl;

    if (!hyperbolic)
    {
        FatalErrorInFunction
            << "The configuration fails the hyperbolicity check. The report is in "
            << "postProcessing/hyperbolicity/report." << nl
            << report.str()
            << exit(FatalError);
    }
}


// ************************************************************************* //
