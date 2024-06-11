//This includes the header file for the OpenFOAM libraries
#include "fvCFD.H"

#include <sstream>

//I0dictionary class is declared in $FOAM_SRC/OpenFOAM/db/I0objects/I0dictionary
//This is our block to create the class ftlDict derived from the IOdictionary class
class ftlDict : public IOdictionary
{
public:
    //This is the constructor of the ftlDict class, it accepts an I0object as a parameter
    ftlDict(const IOobject& ioObj);

    //This is the destructor of the ftlDict class
    ~ftlDict();

    void printTokensInTheDict() const;
};