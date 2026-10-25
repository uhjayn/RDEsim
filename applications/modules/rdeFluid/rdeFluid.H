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

    // oneStepProgress or twoStep. oneStepProgress is the default.
    word reactionModel_;

    volScalarField lambda_;

    // twoStep: induction progress xi (0 to 1) and its Arrhenius-type constants.
    autoPtr<volScalarField> xi_;
    scalar inductionA_;
    scalar inductionTa_;
    scalar releaseA_;
    scalar releaseTa_;
    scalar releaseOrder_;
    scalar inductionMargin_;

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

    // Optional. When enabled, a cell whose kinetic energy would drive T
    // below positivityTMin_ has |U| scaled down so e + 1/2|U|^2 is unchanged.
    // The count is the number of cell-steps, not a silent clamp.
    bool positivityEnabled_;
    scalar positivityTMin_;
    label nPositivity_;
    label nEnergyFloor_;


    tmp<volScalarField> reactionRate() const;

    void inductionEquation();

    void progressEquation();

    void boundPositivity();

    void writePositivityCount() const;

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
