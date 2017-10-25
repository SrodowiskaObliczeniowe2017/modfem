#ifndef MMPH_MESH_H_
#define MMPH_MESH_H_

#include <memory>
#include <vector>
#include <set>
#include <map>

#include "mmh_intf.h"
#include "uth_log.h"

//#include "fpcm/distributed_mesh.hpp"
#include "fpcm/compressed_mesh.hpp"
#include "EntitiesOwnership.hpp"
#include "mmph_ipid.h"


namespace mmpt {
    class Transferer;
    class TransferOrder;
    class TransferResult;
    class TransferWithOwnershipResult;
    class CompositeTransferWithOwnershipResult;
}

using namespace mmpt;



// typedef std::unordered_map<int,mmpt_ipid> mmpt_owners_map;

//using namespace fpcm::DistributedMesh;


/// mmpt_mesh is responsible for realization of mmph interface.
/// It uses:
/// - ddl_metis_parmetis for initial mesh partition and repartition (if needed)
/// - fpcm::DistributedMesh for distributed computation mesh logic,
/// - mmpt::TransferOrder for encapsulation of cross processes mesh changes decisions
/// - mmpt::TransferResult for storing mesh transfers effects
/// - mmpt::Transferer for realization of cross processes mesh operations
/// mmpt_mesh also coordinates and supports the elements above,
/// and applies effective changes on local mesh data structures through mmh interafce.

struct ApplyChanges {
    std::vector<int> vertices_old2newIds,
    element_pos2newIds;

    std::vector<TransferOrder>  transfer_directions;
};


typedef EntitiesOwnership<mmpt_ipid> EntOwnerships;
typedef EntitiesOwnership<mmpt_ipid>::EntOwn EntOwn;
// Local ID.
typedef int LID;
// Global ID.
typedef fpcm::CompressedMesh::GID GID;


class mmpt_mesh{
    static const int N_TYPES = MMC_ALL_N_TYPES+1;
public:

  int mesh_id;

  Transferer* transferer;

  EntOwnerships ownership;


  // Constructor
  mmpt_mesh() ;

  mmpt_mesh(const mmpt_mesh & other) ;

  ~mmpt_mesh();

  void init();

  void clear();

  const ApplyChanges& apply(const TransferOrder &o, const TransferResult &r);
  const ApplyChanges& apply(const TransferOrder &o, const TransferWithOwnershipResult &r);
  const ApplyChanges& apply(const TransferOrder &o, const CompositeTransferWithOwnershipResult &r);

  void findContestedExternalVerticesOfElements(const std::vector<int> & element_ids, std::vector<int> & externalVertices) const;
  void findContestedVerticesAssigmentToProcs(const std::vector<TransferOrder> &orders, // in:
                                             std::map<int, mmpt_ipid> &contested_vrts_glob_ids) const;  // out

  typedef std::set<int> NEIG_SET;
  NEIG_SET   neighbours;

  template<int type>
  void GetOffspring(           std::vector<EntOwn> & parents,
                               std::vector<LID> & offspring,
                               std::vector<GID> & offspring_gids) const ;

  template<int type>
  void SetOffspring(const   std::vector<EntOwn> & parents,
                    const   std::vector<LID> & offspring,
                    const   std::vector<GID> & offspring_gids);


  bool  check() const;


  GID GlobID_Node(const LID Loc_id) const {
      double coords[3];
      mmr_node_coor(mesh_id,Loc_id,coords);
      return glob_ids.encode(coords);
  }

private:
    std::vector<ApplyChanges*>   history;
    fpcm::CompressedMesh::CoordMesh glob_ids;

} ;




#endif //MMPH_MESH_H_
