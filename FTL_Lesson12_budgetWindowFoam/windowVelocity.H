forAll(cellExtents,label)
{
    int cellI = cellExtents(label);
    velocityStacks[cellI].push(mag(U[cellI]));

    if(cellI==cellExtents[0])
    {
        Info << "VelocityStacks agter push = "<< velocityStacks[cellI] <<nl <<endl;
    }
    
    if (velocityStacks[cellI].size() == window)
    {
        scalar meanValue = 0;
        forAllConstIter(FIFOStack<scalar>, velocityStacks[cellI],iter)
        {
            meanValue += iter();
        }
        meanValue /= window;
        velocityStacks[cellI].pop();
        velocityWindow[cellI].append(meanValue);
       
        //Uwindowed[cellI] = velocityWindow[cellI].back(); v2306

        if(cellI == cellExtents[0])
        {
            Info << "velocityStacks after Pop = " << velocityStacks[cellI] << nl << endl;
            Info << "velocityWindow + meanValue = " << velocityWindow[cellI] << nl << endl;
            //Info << "Uwindowed = " << Uwindowed[cellI] << nl << endl; v2306
        }
    }
}