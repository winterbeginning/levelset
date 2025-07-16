for (int corr=0; corr<int(epsilon.value()/deltaTau.value()); corr++)
{
    psi = psi + psi0/Foam::sqrt(Foam::sqr(psi0)+Foam::sqr(epsilon))
        *(double(1)-mag(fvc::grad(psi)*dimChange))*deltaTau;   
    psi.correctBoundaryConditions();                                                 
}
