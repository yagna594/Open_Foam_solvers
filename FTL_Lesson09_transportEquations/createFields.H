// create passive scalar Field, T, and read BCs and the initial conditions
Info << "Reading in field T" <<nl <<endl;
volScalarField T
(   
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
    

Info << "Reading in field U" <<nl <<endl;
volVectorField U
(   
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//Read in the transport properties and get diffusion coefficient for T
Info << "Reading in transportProperties" <<nl <<endl;
IOdictionary transportProperties
(   
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::AUTO_WRITE
    )
);

Info<< "Reading in diffusivity" <<nl << endl;
dimensionedScalar DT (transportProperties.lookup("DT"));

//#include "createPhi.H"
Info << "Reading/creating face flux field phi" <<nl <<endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::interpolate(U) & mesh.Sf()   //interpolation scheme specified in fvSchems [m/s]*m2 == m3/s == flow rate
);

