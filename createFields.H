//Start of sclscvof
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar deltaX  //网格尺度
(
    transportProperties.lookup("deltaX")
);

dimensionedScalar gamma  //控制界面过渡区的厚度
(
    dimensionedScalar(deltaX*0.75)
);

dimensionedScalar epsilon  //光滑带宽
(  //        dimensionedScalar(deltaX*1.5)
    transportProperties.lookup("epsilon")
);

dimensionedScalar deltaTau  //伪时间步长
(
    dimensionedScalar(deltaX*0.1)
);

dimensionedScalar dimChange
(
    dimensionedScalar("dimChange",dimLength, 1.0)
);

dimensionedScalar dimtime
(
    dimensionedScalar("dimtime",dimTime, 1.0)
);

// 创建Level-Set场
volScalarField psi
(
    IOobject
    (
        "psi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
volScalarField psi0
(
    IOobject
    (
        "psi0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    psi
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
    //dimensionedVector("U", dimVelocity, vector(0.1,0,0))
);
