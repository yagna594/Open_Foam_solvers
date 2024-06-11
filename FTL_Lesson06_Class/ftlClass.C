//This line includes the header file ftlClass.H which contains the declaration of the
// ftlclass.
#include "ftlClass.H"

// Constructors
//----------------------------------
ftlClass :: ftlClass()
{
    myInt_= 0;
}
// Destructors
//----------------------------------
ftlClass ::~ ftlClass()
{
    // Cleanup code if necessary
}
// Public Member Functions
//---------------------------------

label ftlClass::doubleStoredValue() const
{
    Info << "Calling ftlClass: doubleStoredValue() " <<endl;
    return myInt_*2;
}

void ftlClass::meshCountFunction(fvMesh& mesh)
{
    Info << "ftlClass recieved the mesh, it has " <<mesh.C().size() << " cells " <<endl;
    myInt_ = mesh.C().size();
}
//Private Member Functions