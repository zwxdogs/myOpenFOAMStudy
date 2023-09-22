/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "myMHDInjection.H"
#include "scalarIOList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::myMHDInjection<CloudType>::myMHDInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    inputFileName_(this->coeffDict().lookup("inputFile")),
    duration_(this->coeffDict().template lookup<scalar>("duration")),
    randomise_(readBool(this->coeffDict().lookup("randomise"))),
    nParticles_(0),
    timestepsPerSecond_
    (
        this->coeffDict().template lookup<scalar>("timestepsPerSecond")
    ),
    injectors_
    (
        IOobject
        (
            inputFileName_,
            owner.db().time().constant(),
            owner.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    currentInjector_(0),
    injectorCells_(0),
    injectorTetFaces_(0),
    injectorTetPts_(0)
{
    duration_ = owner.db().time().userTimeToTime(duration_);

    // Set/cache the injector cells
    injectorCells_.setSize(injectors_.size());
    injectorTetFaces_.setSize(injectors_.size());
    injectorTetPts_.setSize(injectors_.size());

    topoChange();

    // Determine volume of particles to inject
    this->volumeTotal_ = 0.0;
    nParticles_ = this->nParticleFixed_;
    forAll(injectors_, i)
    {
        this->volumeTotal_ += injectors_[i].parcelsNeedInject()*(nParticles_)*
                                ((pi)*pow3(injectors_[i].d())/6)*(timestepsPerSecond_);
    }
    this->volumeTotal_ *= duration_;
}


template<class CloudType>
Foam::myMHDInjection<CloudType>::myMHDInjection
(
    const myMHDInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    inputFileName_(im.inputFileName_),
    duration_(im.duration_),
    randomise_(im.randomise_),
    nParticles_(im.nParticles_),
    timestepsPerSecond_(im.timestepsPerSecond_),
    injectors_(im.injectors_),
    currentInjector_(im.currentInjector_),
    injectorCells_(im.injectorCells_),
    injectorTetFaces_(im.injectorTetFaces_),
    injectorTetPts_(im.injectorTetPts_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::myMHDInjection<CloudType>::~myMHDInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::myMHDInjection<CloudType>::topoChange()
{
    // Set/cache the injector cells
    forAll(injectors_, i)
    {
        this->findCellAtPosition
        (
            injectorCells_[i],
            injectorTetFaces_[i],
            injectorTetPts_[i],
            injectors_[i].x()
        );
    }
}


template<class CloudType>
Foam::scalar Foam::myMHDInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::label Foam::myMHDInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    scalar parcelsToInject = 0;

    if ((time0 >= 0.0) && (time0 < duration_))
    {
        /* return floor(injectorCells_.size()*(time1 - time0)*parcelsPerSecond_); */
        forAll(injectors_, i)
        {
            parcelsToInject += injectors_[i].parcelsNeedInject();
        }
    }
    else
    {
        return 0;
    }

    return parcelsToInject;
}


template<class CloudType>
Foam::scalar Foam::myMHDInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    scalar volume = 0.0;
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        forAll(injectors_, i)
        {
            volume += injectors_[i].parcelsNeedInject()*(nParticles_)*
                        ((pi)*pow3(injectors_[i].d())/6)*(timestepsPerSecond_)*(time1 - time0);
        }
    }

    return volume;
}


template<class CloudType>
void Foam::myMHDInjection<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label nParcels,
    const scalar time,
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
)
{
    label injectorI = 0;
    if (randomise_)
    {
        Random& rnd = this->owner().rndGen();
        injectorI = rnd.sampleAB<label>(0, injectorCells_.size());
    }
    else
    {
        injectorI = int64_t(parcelI)*int64_t(injectors_.size())/nParcels;
    }

    position = injectors_[injectorI].x();
    cellOwner = injectorCells_[injectorI];
    tetFacei = injectorTetFaces_[injectorI];
    tetPti = injectorTetPts_[injectorI];
}


template<class CloudType>
void Foam::myMHDInjection<CloudType>::setProperties
(
    const label parcelI,
    const label nParcels,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    label injectorI = int64_t(parcelI)*int64_t(injectors_.size())/nParcels;

    // set particle velocity
    parcel.U() = injectors_[injectorI].U();

    // set particle diameter
    parcel.d() = injectors_[injectorI].d();

    // set particle density
    parcel.rho() = injectors_[injectorI].rho();

    // currentInjector
    currentInjector_ = injectorI;
}

template<class CloudType>
bool Foam::myMHDInjection<CloudType>::fullyDescribed() const
{
    return true;
}


template<class CloudType>
bool Foam::myMHDInjection<CloudType>::validInjection
(
    const label
)
{
    return true;
}


// ************************************************************************* //
