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

    psi = psi0;

    #include "init.H"
    
    // 时间循环
    while (runTime.run())
    {
        runTime++;
        Info << "Time = " << runTime.timeName() << nl << endl;

        volVectorField gradPsi(fvc::grad(psi));   
        //volVectorField nVecfv(gradPsi/(mag(gradPsi)+scalar(1.0e-10)/dimChange)); 
        //volScalarField U_n(U & nVecfv);
        //volVectorField U_nVecfv(U_n * nVecfv);
        volScalarField U_n(U & gradPsi);

        //surfaceScalarField phi_n = fvc::flux(U_nVecfv);

        fvScalarMatrix psiEqn
        (
            fvm::ddt(psi)
          //+ fvm::div(phi_n, psi)
          + U_n
        );

        

        psi = psi0;
        #include "init.H"
        psiEqn.solve();

        runTime.write();
    }

    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
