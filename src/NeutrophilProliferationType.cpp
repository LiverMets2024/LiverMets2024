#include "NeutrophilProliferativeType.hpp"

NeutrophilProliferativeType::NeutrophilProliferativeType()
    : AbstractCellProliferativeType(13)
{}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(NeutrophilProliferativeType)