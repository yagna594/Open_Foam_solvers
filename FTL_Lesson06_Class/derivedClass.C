//This includes the header file derivedClass.H, which contains the declatarion of the
//ftlDict class.
#include "derivedClass.H"

// Constructors
//--------------------------------------------
//This block defines the constructor for the ftlDict class
ftlDict::ftlDict(const IOobject& ioObj)
:
    IOdictionary(ioObj)
{
//Nothing further to initialise since everything is done using the base class
//I0dictionary
}

// Destructors
//--------------------------------------------
//This block defines the destruct
ftlDict::~ ftlDict ()
{}

// Public Member Functions
//--------------------------------------------
void ftlDict::printTokensInTheDict() const
{
    List<token> characters(this->tokens());

    std::stringstream ss;
    ss << "Tokens in the file: ";

    forAll(characters, i)
        if(characters[i].isWord())
            ss << "\n" << tab <<characters[i].wordToken();

    Info << ss.str().c_str() << endl;

}

// Private Member Functions

//If we had private memebr functions they would go here