/*

Copyright (c) 2005-2020, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "AbstractChemokinePde.hpp"
#include "ApoptoticCellProperty.hpp"
#include "CellLabel.hpp"

template<unsigned DIM>
AbstractChemokinePde<DIM>::AbstractChemokinePde(AbstractCellPopulation<DIM,DIM>& rCellPopulation,
                                                            double duDtCoefficient,
                                                            double diffusionCoefficient,
                                                            double decayRate)
    : AveragedSourceParabolicPde<DIM>(rCellPopulation,duDtCoefficient,diffusionCoefficient,0.0),
    m_decayRate(decayRate)
{
}

template<unsigned DIM>
void AbstractChemokinePde<DIM>::SetupSourceTerms(TetrahedralMesh<DIM,DIM>& rCoarseMesh, std::map<CellPtr, unsigned>* pCellPdeElementMap) // must be called before solve
{
    // Allocate memory
    this->mCellDensityOnCoarseElements.resize(rCoarseMesh.GetNumElements());
    this->mCellDensityOnCoarseElements1.resize(rCoarseMesh.GetNumElements());
    for (unsigned elem_index=0; elem_index<this->mCellDensityOnCoarseElements.size(); elem_index++)
    {
        this->mCellDensityOnCoarseElements[elem_index] = 0.0;
        this->mCellDensityOnCoarseElements1[elem_index] = 0.0;
    }

}

template<unsigned DIM>
double AbstractChemokinePde<DIM>::ComputeDuDtCoefficientFunction(const ChastePoint<DIM>& )
{
    return this->mDuDtCoefficient;
}

template<unsigned DIM>
double AbstractChemokinePde<DIM>::ComputeSourceTerm(const ChastePoint<DIM>& rX, double u, Element<DIM,DIM>* pElement)
{
    assert(!this->mCellDensityOnCoarseElements.empty());
    assert(!this->mCellDensityOnCoarseElements1.empty());
    double source_coefficient = this->mCellDensityOnCoarseElements[pElement->GetIndex()];
    double consumption_coefficient = this->mCellDensityOnCoarseElements1[pElement->GetIndex()];
    // The source term is rho(x)*mSourceCoefficient - mDecayRate*u
    return source_coefficient + consumption_coefficient * u - m_decayRate * u;
}

// LCOV_EXCL_START
template<unsigned DIM>
double AbstractChemokinePde<DIM>::ComputeSourceTermAtNode(const Node<DIM>& rNode, double u)
{
    NEVER_REACHED;
    return 0.0;
}
// LCOV_EXCL_STOP

template<unsigned DIM>
c_matrix<double,DIM,DIM> AbstractChemokinePde<DIM>::ComputeDiffusionTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement)
{
    return this->mDiffusionCoefficient*identity_matrix<double>(DIM);
}

template<unsigned DIM>
double AbstractChemokinePde<DIM>::GetUptakeRateForElement(unsigned elementIndex)
{
    this->mCellDensityOnCoarseElements[elementIndex];
    return this->mCellDensityOnCoarseElements[elementIndex];
}

template <unsigned DIM>
double AbstractChemokinePde<DIM>::GetDecayRate() const
{
    return m_decayRate;
}

// Explicit instantiation
template class AbstractChemokinePde<1>;
template class AbstractChemokinePde<2>;
template class AbstractChemokinePde<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AbstractChemokinePde)
