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

#include "CXCL8Pde.hpp"

#include "MetCellProliferativeType.hpp"
#include "NeutrophilProliferativeType.hpp"
#include "ApoptoticCellProperty.hpp"
#include "CellLabel.hpp"

template<unsigned DIM>
CXCL8Pde<DIM>::CXCL8Pde(AbstractCellPopulation<DIM,DIM>& rCellPopulation, double duDtCoefficient, double diffusionCoefficient,
        double decayRate,
        double sourceCoefficient,
        double consumptionCoefficient) :
    AbstractChemokinePde<DIM>(rCellPopulation,duDtCoefficient,diffusionCoefficient, decayRate),
    m_sourceCoefficient(sourceCoefficient),
    m_consumptionCoefficient(consumptionCoefficient)
{
}

template<unsigned DIM>
void CXCL8Pde<DIM>::SetupSourceTerms(TetrahedralMesh<DIM,DIM>& rCoarseMesh, std::map<CellPtr, unsigned>* pCellPdeElementMap) // must be called before solve
{
    AbstractChemokinePde<DIM>::SetupSourceTerms(rCoarseMesh, pCellPdeElementMap);

    // Loop over cells, find which coarse element it is in, and add 1 to mSourceTermOnCoarseElements[elem_index]
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mrCellPopulation.Begin();
         cell_iter != this->mrCellPopulation.End();
         ++cell_iter)
    {
        unsigned elem_index = 0;
        const ChastePoint<DIM>& r_position_of_cell = this->mrCellPopulation.GetLocationOfCellCentre(*cell_iter);

        if (pCellPdeElementMap != nullptr)
        {
            elem_index = (*pCellPdeElementMap)[*cell_iter];
        }
        else
        {
            elem_index = rCoarseMesh.GetContainingElementIndex(r_position_of_cell);
        }

        if (!cell_iter->template HasCellProperty<ApoptoticCellProperty>())
        {
            // Count Source Cells
            if (cell_iter->GetCellProliferativeType()->template IsType<MetCellProliferativeType>())
            {
                this->mCellDensityOnCoarseElements[elem_index] += 1.0;
            }
            // Count Consumption Cells
            if (cell_iter->GetCellProliferativeType()->template IsType<NeutrophilProliferativeType>())
            {
                this->mCellDensityOnCoarseElements1[elem_index] += 1.0;
            }
        }
    }

    // Then divide each entry of mSourceTermOnCoarseElements by the element's area
    c_matrix<double, DIM, DIM> jacobian;
    double det;
    for (unsigned elem_index=0; elem_index<this->mCellDensityOnCoarseElements.size(); elem_index++)
    {
        rCoarseMesh.GetElement(elem_index)->CalculateJacobian(jacobian, det);
        this->mCellDensityOnCoarseElements[elem_index] /= rCoarseMesh.GetElement(elem_index)->GetVolume(det);
        this->mCellDensityOnCoarseElements1[elem_index] /= rCoarseMesh.GetElement(elem_index)->GetVolume(det);
    }
}


// Explicit instantiation
template class CXCL8Pde<1>;
template class CXCL8Pde<2>;
template class CXCL8Pde<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CXCL8Pde)