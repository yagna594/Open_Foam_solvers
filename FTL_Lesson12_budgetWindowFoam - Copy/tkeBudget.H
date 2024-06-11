if (U.db().foundObject<volVectorField>("UMean"))
{
    const objectRegistry& db = U.db();
    //k = 1/2(u'u'+ v'v'+ w'w')
    //u' = u-u mean
    // volScalarField kRes = 0.5*tr(db.lookup0bject<volSymmTensorField>("UFrime2Mean")); //Instaneous resolved tke
    volVectorField UPrime = (U - db.lookupObject<volVectorField>("UMean")); //U' = U - U_mean
    volScalarField kRes = 0.5*(UPrime & UPrime) ; //Instantaneous resolved tke
    volScalarField kSGS = turbulence->k();

    volScalarField kTot = kRes + kSGS;

    LESResIndex = kRes/max(kTot,kSmall);

    Info << "ESTIMATED LESResIndex min, max, avg:" << min(LESResIndex).value()
        << ", " << max(LESResIndex).value()
        << ", " << LESResIndex.weightedAverage(mesh.V()).value() << endl;
}

else
{
    Warning << "Averaging of velocity U is not active in the condtrolDict." << endl
            << "This is necessary to calculate resolved TKE" << endl
            << "Add field U to the fieldAverage untility in controlDict." << endl;
}