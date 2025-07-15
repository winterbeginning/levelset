for (int corr=0; corr<int(epsilon.value()/deltaTau.value()); corr++)
{
    psi = psi + psi0/Foam::sqrt(Foam::pow(psi0,2)+Foam::pow(epsilon,2))
        *(double(1)-mag(fvc::grad(psi)*dimChange))*deltaTau;   
    psi.correctBoundaryConditions();                                                 
}

// volScalarField S = psi0/Foam::sqrt(Foam::pow(psi0,2)+Foam::pow(epsilon,2));

// // 3.2 重新初始化迭代（半隐式格式）
// for (int corr=0; corr<int(epsilon.value()/deltaTau.value()); corr++)
// {
//     fvScalarMatrix reinitEqn
//     (
//         fvm::ddt(psi) == S*(1.0 - mag(fvc::grad(psi))*dimChange)/dimtime
//     );
//     reinitEqn.solve();
// }