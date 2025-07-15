    // 创建Level-Set场
    volScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("U", dimVelocity, vector::zero)
    );

    const scalar pi = constant::mathematical::pi;
    forAll(mesh.C(), cellI)
    {
        const scalar x = mesh.C()[cellI].x();
        const scalar y = mesh.C()[cellI].y();
        U[cellI].x() = -Foam::pow(Foam::sin(pi*x), 2) * Foam::sin(2*pi*y);
        U[cellI].y() = -Foam::pow(Foam::sin(pi*y), 2) * Foam::sin(2*pi*x);
        phi[cellI] = Foam::sqrt(Foam::pow(x-0.5, 2) + Foam::pow((y-0.75),2))-0.015;
    }
    U.correctBoundaryConditions();

    // 设置通量场
    surfaceScalarField phiFlux
    (
        IOobject("phiFlux", runTime.timeName(), mesh),
        fvc::flux(U)
    );