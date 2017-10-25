#ifndef ENTITIES_OWNERSHIP_HPP_
#define ENTITIES_OWNERSHIP_HPP_

#include <vector>

#include "uth_log.h"
#include "mmh_intf.h"
#include "un_ord_map_selector.hpp"


namespace mmpt {


typedef int (*fStatPtr)(int,int);
typedef int (*fQueryPtr)(int);

const int N_TYPES = MMC_ALL_N_TYPES+1;

static const fStatPtr status_func[N_TYPES] =
{
    mmr_node_status,
    mmr_edge_status,
    mmr_fa_status,
    mmr_el_status,
    NULL
};

static const fQueryPtr last_func[N_TYPES] =
{
    mmr_get_max_node_id,
    mmr_get_max_edge_id,
    mmr_get_max_face_id,
    mmr_get_max_elem_id,
    NULL
};


template<typename TOwnership>
class EntitiesOwnership
{
    int local_owner_id;
    int mesh_id;

public:

    typedef TOwnership EntOwn;


    static const EntOwn OWN_NULL ;
    static const EntOwn OWN_SELF ;


    enum BOUNDARABLE_TYPE{
        INNER_BOUNDARABLE=0, // our ownership, but other processes have them
        OUTER_BOUNDARABLE,   // overlap = foregin ownership
        ANY_BOUNDARABLE,
        LAST_BOUNDARABLE
    };

    static const int UNDEFINED = -1;
    static const int SELF = -2;
    static const int FOREGIN = -3;

    EntitiesOwnership()
    : local_owner_id(-1),mesh_id(-1)
    {

    }

    void setMeshId(const int Mesh_id) {
        mesh_id = Mesh_id;

        boundarables[MMC_NODE].resize(mmr_get_nr_node(Mesh_id));
        boundarables[MMC_EDGE].resize(mmr_get_nr_edge(Mesh_id));
        boundarables[MMC_FACE].resize(mmr_get_nr_face(Mesh_id));
        boundarables[MMC_ELEMENT].resize(mmr_get_nr_elem(Mesh_id));
    }

    void setLocalOwnerID(int owner_proc_id) {
        local_owner_id = owner_proc_id;
    }

    int getLocalOwnerID() const {
        return local_owner_id;
    }

    void Set(EntOwn & gid, const int owner, const int id_at_owner) const {
        gid.owner = owner;
        gid.id_at_owner = id_at_owner;
    }

    void SetOwner(EntOwn & gid, const int owner) const {
        gid.owner = owner;
    }

    void SetId_at_owner(EntOwn & gid, const int id_at_owner) const {
        gid.id_at_owner = id_at_owner;
    }

    int GetOwner(const EntOwn & gid) const {
        return gid.owner;
    }

    int GetId_at_owner(const EntOwn & gid) const {
        return gid.id_at_owner;
    }

    template<int type>
    void    AddOwnership(const int loc_id,const EntOwn& own)  {
        mf_check_debug(loc_id > 0, "Wrong local id(%d)!",loc_id);
        mf_check(own.owner != local_owner_id, "Adding ownership for self!");
        ownerships[type][loc_id] = own;
        own2loc[type][own] = loc_id;
        mf_log_info("Adding ownership of %d as (%d,%d)",
                    loc_id,own.owner,own.id_at_owner);
        setBoundarable<type>(loc_id);
    }

    template<int type>
    int GetOwner(const int loc_id) const {
        int owner = UNDEFINED;
        mf_check_debug(loc_id > 0, "Wrong local id(%d)!",loc_id);
        typename O_MAP::const_iterator it = ownerships[type].find(loc_id);
        if(it != ownerships[type].end()) {
            owner = it->second.owner;
        }
        else {
            if(MMC_ACTIVE == status_func[type](mesh_id,loc_id)) {
                owner = local_owner_id;
            }
        }
        return owner;
    }

    template<int type>
    int GetId_at_owner(const int loc_id) const  {
        int id_at_owner = UNDEFINED;
        mf_check_debug(loc_id > 0, "Wrong local id(%d)!",loc_id);
        typename O_MAP::const_iterator it = ownerships[type].find(loc_id);
        if(it != ownerships[type].end()) {
            id_at_owner = it->second.id_at_owner;
        }
        else {
            if(MMC_ACTIVE == status_func[type](mesh_id,loc_id)) {
                id_at_owner = loc_id;
            }
        }
        return id_at_owner;
    }

    template<int type>
    void GetOwnership(const int loc_id, EntOwn& own) const {
        mf_check_debug(loc_id > 0, "Wrong local id(%d)!",loc_id);
        typename O_MAP::const_iterator it = ownerships[type].find(loc_id);
        if(it == ownerships[type].end() ) {
            own.id_at_owner = loc_id;
            own.owner = local_owner_id;
        }
        else {
            own = it->second;
        }
    }

    template<int type>
    int GetLocIdForOwnership(const EntOwn & own) const {
        int loc_id = UNDEFINED;
        if(own.owner == local_owner_id) {
            if(MMC_ACTIVE == status_func[type](mesh_id,own.id_at_owner)) {
                loc_id = own.id_at_owner;
            }
        }
        else {
            typename G_MAP::const_iterator it = own2loc[type].find(own);
            if(it != own2loc[type].end()) {
                loc_id = it->second;
            }
        }

        mf_check(loc_id <= last_func[type](mesh_id), "Loc id (%d) out-of-bounds(%d)!",loc_id, last_func[type](mesh_id));

        return loc_id;
    }

    template<int type>
    int RemoveID(const int LocID)
    {
        typename O_MAP::iterator it = ownerships[type].find(LocID);
        if(it != ownerships[type].end()) {
            own2loc[type].erase(it->second);
            ownerships[type].erase(it);
        }
        boundarables[type][LocID]=false;
    }

    template<int type>
    int RemoveID(const EntOwn & own)
    {
        typename G_MAP::iterator it = own2loc[type].find(own);
        if(it != own2loc[type].end()) {
            ownerships[type].erase(it->second);
            own2loc[type].erase(it);
        }
    }

    template<int type>
    void setBoundarable(const int Loc_pos)
    {
        if(boundarables[type].size() <= Loc_pos) {
            boundarables[type].resize(Loc_pos+10);
        }
        boundarables[type][Loc_pos] = true;
    }

    template<int type>
    const std::vector<bool>& getBoundarables() const {
        return boundarables[type];
    }


    template<int Ttype>
    void getBoundarablesGlobIDs(std::vector<EntOwn> & Adapt_boundarables,
                                     const BOUNDARABLE_TYPE B_type) const
    {
        std::vector<bool>::const_iterator it = boundarables[Ttype].begin();
        const std::vector<bool>::const_iterator end = boundarables[Ttype].end();

        for(int i=0;it != end;++it, ++i) {
            if(*it) {
                typename O_MAP::const_iterator o_it = ownerships[Ttype].find(i);

                bool is_local = ( o_it == ownerships[Ttype].end() ) ;

                switch(B_type) {
                case INNER_BOUNDARABLE:
                    if(is_local) {
                         Adapt_boundarables.push_back(EntOwn(local_owner_id,i));
                    }
                    break;
                case OUTER_BOUNDARABLE:
                    if(!is_local) {
                        Adapt_boundarables.push_back(o_it->second);
                    }
                    break;
                case ANY_BOUNDARABLE:
                default:
                    mf_fatal_err("Unknown boundrable type");
                    break;
                }
            }
        }
    }

    template<int type>
    int getBoundaryOwnerships(int * first, const int *last,
                       std::vector<const EntOwn*>& boundOwns) const
    {
        assert(first <= last);

        int n_boundarables=0;

        boundOwns.reserve(last-first);

        typename O_MAP::const_iterator it;

        for(;first != last; ++first) {
            if(boundarables[type][*first]) {
                it = ownerships[type].find(*first);
                if(it == ownerships[type].end()) {
                    boundOwns.push_back(& OWN_SELF);
                }
                else {
                    boundOwns.push_back( & it->second );
                }

                ++n_boundarables;
            }
            else {
                boundOwns.push_back(& OWN_NULL);
            }
        }
        return n_boundarables;
    }

    bool check() const
    {
        int nel=0;
        EntOwn tmp;
        while(0!=(nel=mmr_get_next_elem_all(mesh_id,nel))) {
            GetOwnership<MMC_ELEMENT>(nel,tmp);
            if(tmp.owner == local_owner_id) {
                mf_check(tmp.id_at_owner == nel,
                         "Incorrect ownership(%d,%d) for %d !",
                         tmp.owner,tmp.id_at_owner,
                         nel);
            }

            mf_check(GetLocIdForOwnership<MMC_ELEMENT>(tmp) == nel,
                     "Incorrect ownership(%d,%d) for %d !",
                     tmp.owner,tmp.id_at_owner,
                     nel);

        }



        return true;
    }

private:

#ifdef BOOST_NO_CXX11_HDR_UNORDERED_MAP
    typedef boost::unordered_map<int,EntOwn> O_MAP;
    typedef boost::unordered_map<EntOwn,int> G_MAP;
#else
    typedef std::unordered_map<int,EntOwn> O_MAP;
    typedef std::unordered_map<EntOwn,int> G_MAP;
#endif

    O_MAP ownerships[N_TYPES];
    G_MAP   own2loc[N_TYPES];

    std::vector<bool> boundarables[N_TYPES];
};

template<typename T>
const typename EntitiesOwnership<T>::EntOwn EntitiesOwnership<T>::OWN_NULL(EntitiesOwnership<T>::UNDEFINED,EntitiesOwnership<T>::UNDEFINED);
template<typename T>
const typename EntitiesOwnership<T>::EntOwn EntitiesOwnership<T>::OWN_SELF(EntitiesOwnership<T>::SELF,EntitiesOwnership<T>::UNDEFINED);


}

#endif // ENTITIES_OWNERSHIP_HPP_
