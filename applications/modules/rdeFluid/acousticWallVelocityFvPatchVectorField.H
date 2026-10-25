/*---------------------------------------------------------------------------*\
  RDEsim acousticWallVelocity
  SPDX-License-Identifier: GPL-3.0-or-later

  Local wall-normal velocity for a microstructure boundary.

  slip       zero normal velocity, tangential velocity from the interior
  resistive  u_n = (p - pRef) / ZReal
  helmholtz  m xi'' + r xi' + s xi = p',  u_face = porosity * xi'

  Coefficients are the values written by the case generator. This condition
  does not invent a cavity geometry, and it does not enforce the imaginary
  part of a quarter-wave or groove impedance.
\*---------------------------------------------------------------------------*/

#ifndef acousticWallVelocityFvPatchVectorField_H
#define acousticWallVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

namespace Foam
{

class acousticWallVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    word model_;
    word pName_;
    scalar pRef_;
    scalar ZReal_;
    scalar resistiveFloor_;
    scalar inertia_;
    scalar resistance_;
    scalar stiffness_;
    scalar porosity_;
    scalar initialDisplacement_;
    scalar initialVelocity_;
    scalarField xi_;
    scalarField v_;


public:

    TypeName("acousticWallVelocity");


    acousticWallVelocityFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, fvMesh>&,
        const dictionary&
    );

    acousticWallVelocityFvPatchVectorField
    (
        const acousticWallVelocityFvPatchVectorField&,
        const fvPatch&,
        const DimensionedField<vector, fvMesh>&,
        const fieldMapper&
    );

    acousticWallVelocityFvPatchVectorField
    (
        const acousticWallVelocityFvPatchVectorField&,
        const DimensionedField<vector, fvMesh>&
    );

    virtual tmp<fvPatchVectorField> clone
    (
        const DimensionedField<vector, fvMesh>& iF
    ) const
    {
        return tmp<fvPatchVectorField>
        (
            new acousticWallVelocityFvPatchVectorField(*this, iF)
        );
    }

    virtual void map(const fvPatchVectorField&, const fieldMapper&);

    virtual void reset(const fvPatchVectorField&);

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

}

#endif

// ************************************************************************* //
