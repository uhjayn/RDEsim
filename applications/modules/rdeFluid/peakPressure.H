/*---------------------------------------------------------------------------*\
  RDEsim peak-pressure function object
  SPDX-License-Identifier: GPL-3.0-or-later

  Stores the running maximum of p in each cell and the time that maximum
  occurred. These fields are a soot-foil record for this one-step model.
  They are not a measured detonation-cell size.
\*---------------------------------------------------------------------------*/

#ifndef peakPressure_H
#define peakPressure_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "zeroGradientFvPatchFields.H"

namespace Foam
{
namespace functionObjects
{

class peakPressure
:
    public fvMeshFunctionObject
{
    volScalarField pMax_;
    volScalarField tAtPMax_;

    bool sample();


public:

    TypeName("peakPressure");

    peakPressure
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    peakPressure(const peakPressure&) = delete;

    virtual ~peakPressure();

    virtual wordList fields() const;

    virtual bool execute();

    virtual bool write();

    void operator=(const peakPressure&) = delete;
};

}
}

#endif

// ************************************************************************* //
