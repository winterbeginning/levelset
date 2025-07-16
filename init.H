for (int corr=0; corr<int(epsilon.value()/deltaTau.value()); corr++)
{
    psi = psi + psi0/Foam::sqrt(Foam::pow(psi0,2)+Foam::pow(epsilon,2))
        *(double(1)-mag(fvc::grad(psi)*dimChange))*deltaTau;   
    psi.correctBoundaryConditions();                                                 
}
