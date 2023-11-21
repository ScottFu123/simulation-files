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

#include "WenYuDragForce.H"
#include "volFields.H"
#include "fvCFD.H"
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

#include "setRootCaseLists.H"
#include "createTime.H"
#include "createMesh.H"

template<class CloudType>
Foam::scalar Foam::WenYuDragForce<CloudType>::CdRe(const scalar Re) const
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
Foam::WenYuDragForce<CloudType>::WenYuDragForce
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
{}


template<class CloudType>
Foam::WenYuDragForce<CloudType>::WenYuDragForce
(
    const WenYuDragForce<CloudType>& df
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
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::WenYuDragForce<CloudType>::~WenYuDragForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

struct faceINfo
{
  double x;
  double y;
  double z;
  double distanceSqure;
}


bool compare(const faceINfo& face1,const faceINfo& face2)

{

return face1.distanceSqure < face2.distanceSqure;

}



template<class CloudType>
Foam::forceSuSp Foam::WenYuDragForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    
    word findAPatch("elecwall");
    label ID=mesh.boundary().findPatchID(findAPatch);
   
    vector<faceINfo>faceLists;

    for(label faseI=0; faseI < mesh.boundary()[ID].Cf().size;faseI++)
	{
	  if(p.z()>-0.0174)
            {
	       double distanceSqure = (mesh.boundary()[ID].Cf().x()-p.x())*(mesh.boundary()[ID].Cf().x()-p.x()) +
				      (mesh.boundary()[ID].Cf().y()-p.y())*(mesh.boundary()[ID].Cf().y()-p.y()) +
				      (mesh.boundary()[ID].Cf().z()-p.z())*(mesh.boundary()[ID].Cf().z()-p.z());
   
               faceINfo signalface{mesh.boundary()[ID].Cf().x(), mesh.boundary()[ID].Cf().y(), mesh.boundary()[ID].Cf().z(), distanceSqure};
	       faceLists.push_back(signalface);
	    }
	}

    sort(faceLists.begin(),faceLists.end(),compare);

    double fx=fy=fz=0;
    double coff=0.0001;
    for(int i=0; i<5; i++)
	{
           fx+=coff/distanceSqure * (faceLists[i].x-p.x());
           fy+=coff/distanceSqure * (faceLists[i].y-p.y());          
           fz+=coff/distanceSqure * (faceLists[i].z-p.z());
	}
 




    scalar alphac(alphac_[p.cell()]);

    forceSuSp value
    (
        Zero,
        (mass/p.rho())
       *0.75*CdRe(alphac*Re)*muc*pow(alphac, -2.65)/(alphac*sqr(p.d()))
    );

    value.Su() = (fx,fy,fz);

    return value;
}


// ************************************************************************* //
