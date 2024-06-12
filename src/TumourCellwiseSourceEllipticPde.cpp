/*

Copyright (c) 2005-2024, University of Oxford.
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

#include "TumourCellwiseSourceEllipticPde.hpp"
#include "MetCellProliferativeType.hpp"

template<unsigned DIM>
TumourCellwiseSourceEllipticPde<DIM>::TumourCellwiseSourceEllipticPde(AbstractCellPopulation<DIM,DIM>& rCellPopulation, double sourceCoefficient)
    : mrCellPopulation(rCellPopulation),
      mSourceCoefficient(sourceCoefficient)
{
}

template<unsigned DIM>
const AbstractCellPopulation<DIM,DIM>& TumourCellwiseSourceEllipticPde<DIM>::rGetCellPopulation() const
{
    return mrCellPopulation;
}

template<unsigned DIM>
double TumourCellwiseSourceEllipticPde<DIM>::GetCoefficient() const
{
    return mSourceCoefficient;
}

template<unsigned DIM>
double TumourCellwiseSourceEllipticPde<DIM>::ComputeConstantInUSourceTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement)
{
    return 0.0;
}

// LCOV_EXCL_START
template<unsigned DIM>
double TumourCellwiseSourceEllipticPde<DIM>::ComputeLinearInUCoeffInSourceTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement)
{
    NEVER_REACHED;
    return 0.0;
}
// LCOV_EXCL_STOP

template<unsigned DIM>
double TumourCellwiseSourceEllipticPde<DIM>::ComputeLinearInUCoeffInSourceTermAtNode(const Node<DIM>& rNode)
{
    double source_coefficient = 0.0;
    //this needs to be edited to be more sensible - currently uses the same/source sink rate
    if (mrCellPopulation.IsPdeNodeAssociatedWithNonApoptoticCell(rNode.GetIndex()))
    {
        if (mrCellPopulation.GetCellUsingLocationIndex(rNode.GetIndex())-> template HasCellProperty<MetCellProliferativeType>())
        {
            source_coefficient = 5*mSourceCoefficient;
        }
        else
        {
            source_coefficient = -mSourceCoefficient;

        }
    }

    return source_coefficient;
}

template<unsigned DIM>
c_matrix<double,DIM,DIM> TumourCellwiseSourceEllipticPde<DIM>::ComputeDiffusionTerm(const ChastePoint<DIM>& rX)
{
    return identity_matrix<double>(DIM);
}

// Explicit instantiation
template class TumourCellwiseSourceEllipticPde<1>;
template class TumourCellwiseSourceEllipticPde<2>;
template class TumourCellwiseSourceEllipticPde<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TumourCellwiseSourceEllipticPde)
