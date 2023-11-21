/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
#include <vector>
#include "ExtraWenYuDragForce.H"
#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::ExtraWenYuDragForce<CloudType>::CdRe(const scalar Re) const
{
    if (Re > 1000.0)
    {
        return 0.44*Re;
    }
    else
    {
        return 24.0*(1.0 + 0.15*pow(Re, 0.687));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class CloudType>
Foam::ExtraWenYuDragForce<CloudType>::ExtraWenYuDragForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    alphac_
    (
        this->mesh().template lookupObject<volScalarField>
        (
            this->coeffs().lookup("alphac")
        )
    )
{
    //Mesh = &mesh;
}


template<class CloudType>
Foam::ExtraWenYuDragForce<CloudType>::ExtraWenYuDragForce
(
    const ExtraWenYuDragForce<CloudType>& df
)
:
    ParticleForce<CloudType>(df),
    alphac_
    (
        this->mesh().template lookupObject<volScalarField>
        (
            this->coeffs().lookup("alphac")
        )
    )
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ExtraWenYuDragForce<CloudType>::~ExtraWenYuDragForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
struct faceInfo
{
   Foam::scalar x;
   Foam::scalar y;
   Foam::scalar z;
   Foam::scalar distanceSqure;
};


inline bool cmp(faceInfo face1, faceInfo face2)
{
   return face1.distanceSqure < face2.distanceSqure;
}


template<class CloudType>
Foam::forceSuSp Foam::ExtraWenYuDragForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    word findAPatch("rock");
    label ID=this->mesh().boundaryMesh().findPatchID(findAPatch);

    std::vector<faceInfo>faceLists;
    
    int Size = 15;
    faceLists.resize(Size);
    
    for(int m=0; m<Size; m++)
    {
          scalar distanceSqure = 
                              ((this->mesh().boundary()[ID].Cf()[m].component(0)-p.position().component(0))
			      *(this->mesh().boundary()[ID].Cf()[m].component(0)-p.position().component(0)) +
                               (this->mesh().boundary()[ID].Cf()[m].component(1)-p.position().component(1))
                              *(this->mesh().boundary()[ID].Cf()[m].component(1)-p.position().component(1)) +
                               (this->mesh().boundary()[ID].Cf()[m].component(2)-p.position().component(2))
                              *(this->mesh().boundary()[ID].Cf()[m].component(2)-p.position().component(2)));

     	faceInfo signalface = { this->mesh().boundary()[ID].Cf()[m].component(0)-p.position().component(0), 
                            	this->mesh().boundary()[ID].Cf()[m].component(1)-p.position().component(1), 
                            	this->mesh().boundary()[ID].Cf()[m].component(2)-p.position().component(2), distanceSqure};
     	faceLists[m] = signalface;
    }
    sort(faceLists.begin(),faceLists.end(),cmp);

    scalar min=0.0049;
    scalar minx=0;
    scalar miny=0;
    scalar minz=0;
    scalar U=sqrt(p.U().component(0)*p.U().component(0)+p.U().component(1)*p.U().component(1)+p.U().component(2)*p.U().component(2));
    if( U>1e-6 && p.position().component(2) > -0.017)
    {
        
        sort(faceLists.begin(),faceLists.end(),cmp);

  	for(label faseI = 0; faseI < this->mesh().boundary()[ID].size(); faseI++)
        {
		
          if(abs((this->mesh().boundary()[ID].Cf()[faseI].component(0)-p.position().component(0)))>0.008||
	     abs((this->mesh().boundary()[ID].Cf()[faseI].component(1)-p.position().component(1)))>0.008||
	     abs((this->mesh().boundary()[ID].Cf()[faseI].component(2)-p.position().component(2)))>0.008)
          {continue;}

          scalar distanceSqure = 
                              ((this->mesh().boundary()[ID].Cf()[faseI].component(0)-p.position().component(0))
			      *(this->mesh().boundary()[ID].Cf()[faseI].component(0)-p.position().component(0)) +
                               (this->mesh().boundary()[ID].Cf()[faseI].component(1)-p.position().component(1))
                              *(this->mesh().boundary()[ID].Cf()[faseI].component(1)-p.position().component(1)) +
                               (this->mesh().boundary()[ID].Cf()[faseI].component(2)-p.position().component(2))
                              *(this->mesh().boundary()[ID].Cf()[faseI].component(2)-p.position().component(2)));
    
          if(distanceSqure<faceLists[Size-1].distanceSqure)
          {
     		faceInfo signalface = { this->mesh().boundary()[ID].Cf()[faseI].component(0)-p.position().component(0), 
                            	        this->mesh().boundary()[ID].Cf()[faseI].component(1)-p.position().component(1), 
                            		this->mesh().boundary()[ID].Cf()[faseI].component(2)-p.position().component(2), distanceSqure};
		faceLists[Size-1] = signalface;
          }
        }

    	scalar fx=0;
    	scalar fy=0;
    	scalar fz=0;

	scalar particleQ=1.6e-5;
	scalar E;


	for(int i=0; i<Size; i++)
	{
		float distance = sqrt(faceLists[i].distanceSqure);
	        E= 5.0*(-0.000001* distance +0.00000007);
	        if(E<0){E=0;}
		fx+=E * particleQ* (faceLists[i].x)/distance;
    		fy+=E * particleQ* (faceLists[i].y)/distance;          
    		fz+=E * particleQ* (faceLists[i].z)/distance;
	}	
	
	scalar alphac(alphac_[p.cell()]);

   	 forceSuSp value
    	(
        	Zero,
        	(mass/p.rho())
       		*0.75*CdRe(alphac*Re)*muc*pow(alphac, -2.65)/(alphac*sqr(p.d()))
    	);

   	 Foam::vector elecforce(fx,fy,fz);
    	value.Su() = elecforce;
    	return value;
     }
     else
     {
	scalar alphac(alphac_[p.cell()]);

   	 forceSuSp value
    	(
        	Zero,
        	(mass/p.rho())
       		*0.75*CdRe(alphac*Re)*muc*pow(alphac, -2.65)/(alphac*sqr(p.d()))
    	);
    	return value;
     }
}


// ************************************************************************* //
