IOobject dictIO
(
    "windowProperties",
    mesh.time().constant(),
    mesh,
    IOobject :: MUST_READ
);
dictionary windowProperties;
windowProperties = IOdictionary(dictIO);

vector extents_min = windowProperties.lookupOrDefault<vector>("extents_min", vector(0,0,0));
vector extents_max = windowProperties. lookupOrDefault<vector>("extents_max", vector(100,100,100));

Info << "Extents min = " << extents_min << " and max = " << extents_max << "\n" << endl;

DynamicList<scalar> cellExtents;

forAll(mesh.C(),cellI)
{

    if(mesh.C()[cellI][0]>extents_min[0] && mesh.C()[cellI][0]<extents_max[0])
    {
        if(mesh.C()[cellI][1]>extents_min[1] && mesh.C()[cellI][1]<extents_max[1])
        {
            if(mesh.C()[cellI][2]>extents_min[2] && mesh.C()[cellI][2]<extents_max[2])
            {
                cellExtents.append(cellI);
            }
        }
    }
}

Info << "Total Number of Cells = " << cellExtents.size() << "/" << mesh.C().size() << "\n"  <<endl;

scalar windowPeriod = windowProperties.lookupOrDefault<scalar>("window_period", 10);
scalar deltaT = runTime.deltaT().value();
scalar window = int(windowPeriod/deltaT);

Info << "Window Period = " << windowPeriod << "s" << nl << endl;
Info << "Number of entries (Window/deltaT) = " << window << nl << endl;



