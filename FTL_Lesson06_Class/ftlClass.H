//This include necessary Openfoam header file
#include "fvCFD.H"

//This block defines a class named ftlClass
class ftlClass
{
// private member
private:
    label myInt_;

public:
    //constructor of the class, this is called when an object is created
    ftlClass();
    
    //Destructor of the class this is called when an object of the class is destroyed   
    ~ftlClass();

    inline label get() const {return myInt_;};

    inline void set(label newInt) {myInt_ = newInt;};

    label doubleStoredValue() const;

    void meshCountFunction(fvMesh& mesh);
};