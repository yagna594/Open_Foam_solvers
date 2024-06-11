// Info<< "Reading dynamic Mesh file" <<nl <<endl;
IOdictionary dynamicMesh
(
    IOobject
    (
        "dynamicMeshDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);