#ifndef ES_HPP
#define ES_HPP

/// \file es.hpp
/// \brief es - Entity System for finite element mesh.
/// \author kazimierz.michalik@agh.edu.pl


#include <cinttypes>
#include <omp.h>
#include <vector>

#include "dbg.h"
#include "compressed_mesh.hpp"

namespace fpcm
{

template< typename F, template<class T, class = std::allocator<T> > class container_type = std::vector >
int container_in_memory_size (const container_type<F> & c) {
    return sizeof(c) + c.size()*sizeof(typename container_type<F>::value_type);
}

//////////////////
///
/// aka Entity System Mesh Implementation with Components, Assemblages, Processors, Entities.
/// Combining mmh and mmph at once... it will be cool ;)
///
///
namespace ES {

typedef int ID;

typedef int GUID;

typedef void* pType;

typedef unsigned char ComponentFlag;

typedef unsigned char Assemblage;

// Components
typedef int Materialable;
//struct Materialable {
//    static const int ComponentNo = 0;
//    int material_id;
//};

struct SubTypeable {
    static const int ComponentNo = 1;
    pType type;
};

struct Hierarchical {
    static const int ComponentNo = 2;
    ID parent;
};

struct Adaptable {
    static const int ComponentNo = 3;
    int8_t n_sons;
    ID first_son;
};
template<int8_t Tmax_vert>
struct Verticable {
    static const int ComponentNo = 4;
    ID v[Tmax_vert];
    //ID* vts_offset;
};

struct Neigborable{
    static const int ComponentNo = 5;
    ID neigs[2];
};

struct CrossLinkable {
    static const int ComponentNo = 6;
    int other_subdomain;
    int elem_on_other_side;
};

//struct Patterable {
//    static const int ComponentNo = 7;
//    int base;
//    int pattern;
//};

template<typename TComp>
ComponentFlag   flagOf() {
    return 1<<TComp::ComponentNo;
}


//// storages
//struct Storages{
//    std::vector<Verticable> verticables;
//    std::vector<Neigborable> neigborables;
//    std::vector<Adaptable>  adaptables;
//};

//using CompressedMesh::PTID;
//struct Mesh {
//    std::vector<PTID>   points;
//    std::vector<GUID>   edges;
//    std::vector<GUID>   faces;
//    std::vector<GUID>   elements;
//};


// Assemblages
//static const Assemblage
//VertexAssemblage = flagOf<Coordinable>() & flagOf<Hierarchical>(), // VertexAssemblage
//EdgeAssemblage = flagOf<Hierarchical>() & flagOf<Adaptable>() & flagOf<Verticable>(),// EdgeAssemblage
//FaceAssemblage = flagOf<Hierarchical>() & flagOf<Adaptable>() & flagOf<Verticable>() & flagOf<SubTypeable>() & flagOf<Neigborable>(), // Face assemblage
//ElementAssemblage = flagOf<Hierarchical>() & flagOf<Adaptable>() & flagOf<Verticable>() & flagOf<SubTypeable>(); // Elements assemblage

//int PatterableProcessor(const ID *toCategorize[], const int n_entities, const int max_length, Patterable* categorized) {
//    for(int i=0; i < n_entities; ++i) {
//        categorized[i].base = std::min_element(toCategorize[i],toCategorize[i]+max_length);
//        categorized[i].pattern = GetPattern(toCategorize[i]+1,max_length-1);
//    }
//}

//inline int GetPattern(const ID* sequence,int length) {
//    int patternID=0;

//    assert(length > 0);

//    while(--length) {
//    patterns[sequence]
//    }

//    return patternID;
//}

//template<typename TComponent>
//TComponent& getAs(GUID id) {
//    storages[TComponent::ComponentNo][GUID];
//}

//GUID    newGUID(const size_t count) {
//    static GUID lastGUID = 1;
//    const GUID ret;
//#pragma omp atomic capture
//    ret=(lastGUID+=count);

//    return (ret-count);
//}

//int CreationProcessor(const GUID* toCreateEntities, const size_t count, const Assemblage& entitiesAssemblage) {
//    GUID* toCreatePtr = toCreateEntities;
//    const GUID* endPtr = toCreateEntities+count;

//    for(;toCreatePtr != endPtr; ++toCreatePtr) {

//    }
//}

//struct AdaptationScheme {
//    size_t n_entities;
//    Adaptable* toAdaptStorage;
//    Assemblage newEntityAssemblage;

//    int newPerAdaptable=8;
//};



//int AdaptationProcessor(const int count)
//{
//    // toAdapt -> undivided + adopted
//    Hierarchical* toAdapt;

//    const Hierarchical* hPtr = undivided.add(count*8);
//    const Adaptable* aPtr = adopted.add(count);

//#pragma omp parallel for
//    for(int a=0;a < n;++a) {
//        hPtr[a].parent = GUID(toAdapt[a]);
//        aPtr[a].first_son = GUID(hPtr[a]);
//        aPtr[a].n_sons =
//    }
//}

//// Processors
//int AdaptationProcessor2(const GUID* toAdaptEntities, GUID* adoptedEntities, GUID* createdEntities, const AdaptationScheme & scheme)
//{
//    const int limit = scheme.n_entities*scheme.newPerAdaptable;
//    GUID baseGUID = newGUID(limit);

//    GUID*   toAdaptPtr = toAdaptEntities;

//    for(int a=0; a < limit; ++a) {
//        getAs<Adaptable>(toAdaptEntities[a]).first_son = baseGUID;
//        getAs<Adaptable>(toAdaptEntities[a]).n_sons = scheme.newPerAdaptable;
//    }
//    return 0;
//}


//// Storages



}
}


#endif // ES_HPP
