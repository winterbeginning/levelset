/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Application
    icoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "fvModels.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    const scalar pi = constant::mathematical::pi;
    forAll(mesh.C(), cellI)
    {
        const scalar x = mesh.C()[cellI].x();
        const scalar y = mesh.C()[cellI].y();
        U[cellI].x() = -Foam::pow(Foam::sin(pi*x), 2) * Foam::sin(2*pi*y);
        U[cellI].y() = -Foam::pow(Foam::sin(pi*y), 2) * Foam::sin(2*pi*x);
        scalar distance = Foam::sqrt(Foam::pow(x-0.5, 2) + Foam::pow((y-0.75),2))-0.15;
        psi[cellI] = distance;
    }

    // psi0 = psi;
    // #include "init.H"
    
    // 时间循环
    while (runTime.run())
    {
        runTime++;
        Info << "Time = " << runTime.timeName() << nl << endl;

        #include "veo.H"

        //volVectorField gradPsi(fvc::grad(psi));   
        //volScalarField Vn(U & gradPsi); 

        surfaceScalarField phiU = fvc::interpolate(U) & mesh.Sf();

        fvScalarMatrix psiEqn
        (
            fvm::ddt(psi)
          //+ Vn * mag(gradPsi)  
          + fvm::div(phiU, psi, "div(phiU,psi)")
          - fvm::Sp(fvc::div(phiU),psi)
        );

        psiEqn.solve();

        // psi0 = psi;
        // #include "init.H"

        runTime.write();
    }

    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
