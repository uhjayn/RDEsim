/*---------------------------------------------------------------------------*\
  RDEsim rdeFluid
  SPDX-License-Identifier: GPL-3.0-or-later

  Reactive Euler / shockFluid module for RDEsim. It links OpenFOAM, so the
  combined work is GPL-3.0-or-later. See the LICENSE file at the repository
  root. This is not a copy of the shockFluid sources.
\*---------------------------------------------------------------------------*/

#ifndef rdeFluid_H
#define rdeFluid_H

#include "shockFluid.H"

namespace Foam
{
namespace solvers
{

class rdeFluid
:
    public shockFluid
{
    // constant/rdeProperties. Every chemistry coefficient is read from there.
    IOdictionary properties_;

    volScalarField lambda_;

    dimensionedScalar A_;
    scalar beta_;
    dimensionedScalar Ta_;
    dimensionedScalar TRef_;
    dimensionedScalar TFloor_;
    dimensionedScalar Tignition_;
    dimensionedScalar lambdaMin_;
    dimensionedScalar lambdaMax_;
    dimensionedScalar q_;
    word sourceIntegration_;
    word flowModel_;
    word equationOfState_;
    word ddtScheme_;


    tmp<volScalarField> reactionRate() const;

    void progressEquation();

    void checkHyperbolicity() const;

    void writeObjectives() const;


public:

    TypeName("rdeFluid");

    rdeFluid(fvMesh& mesh);

    rdeFluid(const rdeFluid&) = delete;

    virtual ~rdeFluid();

    virtual void prePredictor();

    virtual void thermophysicalPredictor();

    virtual void postSolve();

    void operator=(const rdeFluid&) = delete;
};

}
}

#endif

// ************************************************************************* //
