#include "FibroblastProliferativeType.hpp"

FibroblastProliferativeType::FibroblastProliferativeType()
    : AbstractCellProliferativeType(14)
{}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(FibroblastProliferativeType)