#ifndef NEUTROPHILPROLIFERATIVETYPE_HPP_
#define NEUTROPHILPROLIFERATIVETYPE_HPP_

#include "AbstractCellProliferativeType.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * Subclass of AbstractCellProliferativeType defining a T cell.
 */
class NeutrophilProliferativeType : public AbstractCellProliferativeType
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell proliferative type.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellProliferativeType>(*this);
    }

public:
    /**
     * Constructor.
     */
    NeutrophilProliferativeType();
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(NeutrophilProliferativeType)

#endif /*NEUTROPHILPROLIFERATIVETYPE_HPP_*/