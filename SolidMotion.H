const labelList& owner = mesh.owner();
const labelList& neighbour = mesh.neighbour();
const volVectorField& C = mesh.C();
const surfaceVectorField& Sf = mesh.Sf(); // Face area vectors

dimensionedScalar dt = runTime.deltaT();

volScalarField& alpha2 = mesh.lookupObjectRef<volScalarField>("alpha2");
volScalarField& alpha1 = mesh.lookupObjectRef<volScalarField>("alpha1");
volScalarField& T = mesh.lookupObjectRef<volScalarField>("T");


// Map: Key = targetCell index, Value = List of (direction vector, face area)
Map<DynamicList<Tuple2<vector, scalar>>> targetData;


forAll(mesh.cells(), celli)
{
    // Check if current cell is a "source" (alpha1 < 1)
    if (alpha1[celli] < 1.0)
    {
        const labelList& cFaces = mesh.cells()[celli];

        forAll(cFaces, i)
        {
            label facei = cFaces[i];

            // Only process internal faces
            if (mesh.isInternalFace(facei))
            {
                label own = owner[facei];
                label nei = neighbour[facei];
                
                // Determine which cell is the 'other' one
                label otherCell = (own == celli) ? nei : own;

                // Check if the neighbor is a "target" (alpha1 == 1)
                if (alpha1[otherCell] == 1.0)
                {
                    // Sf always points owner -> neighbour. 
                    // If current cell is owner, Sf is already outward.
                    vector outwardSf = (own == celli) ? Sf[facei] : -Sf[facei];
                    scalar S = mag(outwardSf);

                    targetData(celli).append(Tuple2<vector, scalar>(outwardSf, S));
                }
            }
        }
    }
}


// Accessing the data
forAllConstIters(targetData, iter)
{
    label targetCellIdx = iter.key();
    const DynamicList<Tuple2<vector, scalar>>& entries = iter.val();
        // Get the volume of the current target cell
    scalar V_target = mesh.V()[targetCellIdx];
	scalar h;
	h=Cp1.value()*T[targetCellIdx]+(alpha2[targetCellIdx]*Hf.value());
	h= (h*V_target);
	Info<<iter.key()<<" "<<h<<endl;
    forAll(entries, i)
    {
        vector d = entries[i].first();
        scalar area = entries[i].second();
        
		// Access Us at the specific cell index 'own' or 'nei' 
        // depending on which cell is the source
        h += (dt.value() *  (Us[targetCellIdx] & d) * ((Cp1.value() * Tmelt.value()) + Hf.value()));
		Info<<(Us[targetCellIdx] & d)<<" "<<h<<endl;

        // Your logic here, e.g., flux calculation:
        // scalar heatFlux = ... * area;
    }

	
	    Info<<alpha2[targetCellIdx]<<endl;
	    T[targetCellIdx] = Tmelt.value();
		alpha2[targetCellIdx] = ((h -V_target*(Cp1.value()*Tmelt.value()))/(Hf.value()*V_target));
		
		if (alpha1[targetCellIdx]<0)
		{
		   alpha1[targetCellIdx]=0.0;	
			
		}
		if (alpha1[targetCellIdx]>1)
		{
		  alpha1[targetCellIdx]=1.0	;
			
		}
		alpha1[targetCellIdx]=scalar(1.0)-alpha2[targetCellIdx];
	
	if(alpha2[targetCellIdx]>0)
	{
	Info<<alpha1[targetCellIdx]<<" " <<h<<" "<<targetCellIdx<<endl;
	}
}
    //Info<< "Cell " << targetCellIdx << " has " << vecs.size() 
        //<< " vectors from alpha1>0 neighbors." << endl;


	// need to consider where the area normal vector comes from,consider velocity dot product,find common area between cells.-done
	//calculate alpha1 and alpha2 again.-done
	
	//may need to calculate the extra enthalpy crossing bounds and limit it cell be cell.-to do
	//contribution may be considered by weighting with respect to dot product value / total velocity dot product from negihbour cells-recursive function.-to do.

//  do a block test.-in progress
	

