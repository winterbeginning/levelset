const scalar t = runTime.value();
forAll(mesh.C(), cellI)
{
    const scalar x = mesh.C()[cellI].x();
    const scalar y = mesh.C()[cellI].y();
    U[cellI].x() = 2. * Foam::cos(pi * t /8.) * (-Foam::sin(pi * x) * Foam::sin(pi * x)) * Foam::sin(pi * y) * Foam::cos(pi * y);
    U[cellI].y() = 2. * Foam::cos(pi * t /8.) * (Foam::sin(pi * y) * Foam::sin(pi * y)) * Foam::sin(pi * x) * Foam::cos(pi * x);
}