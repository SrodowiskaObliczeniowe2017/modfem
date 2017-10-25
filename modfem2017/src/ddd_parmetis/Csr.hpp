#ifndef _CSR_HPP_
#define _CSR_HPP_


#include <vector>
#include <iterator>
#include <fstream>
#include <utility>
#include <iostream>

#include "mmh_intf.h"
#include "uth_log.h" //#include "dbg.h"

using namespace std::rel_ops;

template<typename T>
class CSR
{
public:
    typedef T Tind;
    typedef Tind* pTind;
    typedef const Tind* cpTind;

    static const int n_max_neig=MMC_MAXELFAC;

    friend class CsrNode;
    friend class CsrNodeInternal;


    //////////////////////////////////////////////////////////////
    typename std::vector<Tind>::iterator    write(std::vector<Tind> & v) const {
//        mf_log_info("Dumping write to file!!!");
//        std::ofstream file("write.txt");
//        std::ostream_iterator<int> out_it (file,"\n");
//        std::copy(xadj_.begin(),xadj_.end(),out_it);
//        std::copy(adjncy_.begin(),adjncy_.end(),out_it);
//        std::copy(vwgt_.begin(),vwgt_.end(),out_it);
//        std::copy(adjwgt_.begin(),adjwgt_.end(),out_it);
//        std::copy(el_loc_id_.begin(),el_loc_id_.end(),out_it);
//        std::copy(part_.begin(),part_.end(),out_it);
//        file.flush();

        v.reserve(size());
        v.push_back(xadj_.size());
        v.insert(v.end(),xadj_.begin(),xadj_.end());
        v.push_back(adjncy_.size());
        v.insert(v.end(),adjncy_.begin(),adjncy_.end());
        v.push_back(vwgt_.size());
        v.insert(v.end(),vwgt_.begin(),vwgt_.end());
        v.push_back(adjwgt_.size());
        v.insert(v.end(),adjwgt_.begin(),adjwgt_.end());
        v.push_back(el_loc_id_.size());
        v.insert(v.end(),el_loc_id_.begin(),el_loc_id_.end());
        v.push_back(part_.size());
        v.insert(v.end(),part_.begin(),part_.end());
        v.push_back(adj_neig_no_.size());
        v.insert(v.end(),adj_neig_no_.begin(),adj_neig_no_.end());

//        mf_log_info("Dumping write_v to file!!!");
//        std::ofstream file3("write_v.txt");
//        for(int i=0; i < v.size(); ++i) {
//            file3 << v[i] << "\n";
//        }
//        file3.flush();

        return v.end();
    }
//    std::ostream& write(std::ostream & os) const {
//        os << size() << " \n";
//        for(const auto & el: xadj_) { os << el<< " "; } os << "\n";
//        for(const auto & el: adjncy_) { os << el<< " "; } os << "\n";
//        for(const auto & el: vwgt_) { os << el<< " "; } os << "\n";
//        for(const auto & el: adjwgt_) { os << el<< " "; } os << "\n";
//        for(const auto & el: el_loc_id_) { os << el<< " "; } os << "\n";
//        for(const auto & el: part_) { os << el<< " "; } os << "\n";
//        for(const auto & el: adj_neig_no_) { os << el<< " "; } os << "\n\n";
////        os.write(reinterpret_cast<const char*>(xadj()),xadj_.size()*sizeof(T));
////        os.write(reinterpret_cast<const char*>(adjncy()),adjncy_.size()*sizeof(T));
////        os.write(reinterpret_cast<const char*>(vwgt()),vwgt_.size()*sizeof(T));
////        os.write(reinterpret_cast<const char*>(adjwgt()),adjwgt_.size()*sizeof(T));
////        os.write(reinterpret_cast<const char*>(el_id()),el_loc_id_.size()*sizeof(T));
////        os.write(reinterpret_cast<const char*>(part()),part_.size()*sizeof(T));
//        return os;
//    }

    typename std::vector<T>::const_iterator    read(const std::vector<T> & v) {
//        mf_log_info("Dumping read_v to file!!!");
//        std::ofstream file3("read_v.txt");
//        for(int i=0; i < v.size(); ++i) {
//            file3 << v[i] << "\n";
//        }
//        file3.flush();

        typename std::vector<T>::const_iterator it(v.begin());
        xadj_.resize(*it); ++it;
        std::copy(it,it+xadj_.size(),xadj_.begin());
        it+=+xadj_.size();

        adjncy_.resize(*it); ++it;
        std::copy(it,it+adjncy_.size(),adjncy_.begin());
        it+=+adjncy_.size();

        vwgt_.resize(*it); ++it;
        std::copy(it,it+vwgt_.size(),vwgt_.begin());
        it+=+vwgt_.size();

        adjwgt_.resize(*it); ++it;
        std::copy(it,it+adjwgt_.size(),adjwgt_.begin());
        it+=+adjwgt_.size();

        el_loc_id_.resize(*it); ++it;
        std::copy(it,it+el_loc_id_.size(),el_loc_id_.begin());
        it+=+el_loc_id_.size();

        part_.resize(*it); ++it;
        std::copy(it,it+part_.size(),part_.begin());
        it+=+part_.size();

        adj_neig_no_.resize(*it); ++it;
        std::copy(it,it+adj_neig_no_.size(),adj_neig_no_.begin());
        it+=+adj_neig_no_.size();


//        mf_log_info("Dumping read to file!!!");
//        std::ofstream file2("read.txt");
//        std::ostream_iterator<int> out_it2 (file2,"\n");
//        std::copy(xadj_.begin(),xadj_.end(),out_it2);
//        std::copy(adjncy_.begin(),adjncy_.end(),out_it2);
//        std::copy(vwgt_.begin(),vwgt_.end(),out_it2);
//        std::copy(adjwgt_.begin(),adjwgt_.end(),out_it2);
//        std::copy(el_loc_id_.begin(),el_loc_id_.end(),out_it2);
//        std::copy(part_.begin(),part_.end(),out_it2);
//        file2.flush();

        return it;
    }

    std::istream& read(std::istream & is) {
        T size=0;
        is >> size;
        resize(size);
        std::istream_iterator<T> iit(is);

        std::copy(iit,iit+xadj_.size(),xadj());
        std::copy(iit,iit+adjncy_.size(),adjncy());
        std::copy(iit,iit+vwgt_.size(),vwgt());
        std::copy(iit,iit+adjwgt_.size(), adjwgt());
        std::copy(iit,iit+el_loc_id_.size(),el_id());
        std::copy(iit,iit+part_.size(), part());
        std::copy(iit,iit+adj_neig_no_.size(),adj_neig_no());
//        is.read(reinterpret_cast<char*>(xadj()),xadj_.size()*sizeof(T));
//        is.read(reinterpret_cast<char*>(adjncy()),adjncy_.size()*sizeof(T));
//        is.read(reinterpret_cast<char*>(vwgt()),vwgt_.size()*sizeof(T));
//        is.read(reinterpret_cast<char*>(adjwgt()),adjwgt_.size()*sizeof(T));
//        is.read(reinterpret_cast<char*>(el_id()),el_loc_id_.size()*sizeof(T));
//        is.read(reinterpret_cast<char*>(part()),part_.size()*sizeof(T));
        return is;
    }

    friend std::ostream& operator<<(std::ostream & os,  CSR<T> const& csr) {
        return csr.write(os);
    }

    friend std::istream& operator>>(std::istream & is, CSR<T> & csr) {
        return csr.read(is);
    }
//////////////////////////////////////////////////////////////
    class CsrNode
    {
    public:

        friend class iterator;
        friend class const_iterator;

        CsrNode() : el_id(0), vwgt(0), n_neighs(0) {
            std::fill(neighs,neighs+n_max_neig,0);
            std::fill(neig_no,neighs+n_max_neig,0);
        }

        pTind end_neighs() const {return neighs+n_neighs;}

        bool    operator==(const CsrNode & other) const { return el_id == other.el_id; }

        Tind el_id,vwgt,n_neighs;
        Tind neighs[n_max_neig], neig_no[n_max_neig];
    };
//////////////////////////////////////////////////////////////
    class CsrNodeInternal
    {
        friend class iterator;
    public:
        CsrNodeInternal(CSR *container = NULL)
            : pxadj(NULL),padjncy(NULL),pvwgt(NULL),
              padjwgt(NULL),pel_id(NULL),pos_(0),padj_neig_no(NULL) {
            if(container != NULL) {
                pxadj =  container->xadj();
                padjncy =  container->adjncy();
                pvwgt =  container->vwgt();
                padjwgt =  container->adjwgt();
                pel_id =  container->el_id();
                ppart =  container->part();
                padj_neig_no = container->adj_neig_no();
            }
        }

        CsrNodeInternal(const CsrNodeInternal &x) {
            pxadj = x.pxadj;
            pvwgt = x.pvwgt;
            pel_id = x.pel_id;
            ppart = x.ppart;
            pos_ = x.pos_;
            padjncy = x.padjncy;
            padjwgt = x.padjwgt;
            padj_neig_no = x.padj_neig_no;
        }

        CsrNodeInternal( CsrNodeInternal &x) {
            pxadj = x.pxadj;
            pvwgt = x.pvwgt;
            pel_id = x.pel_id;
            ppart = x.ppart;
            pos_ = x.pos_;
            padjncy = x.padjncy;
            padjwgt = x.padjwgt;
            padj_neig_no = x.padj_neig_no;
//            x.pxadj = NULL;
//            x.pvwgt = NULL;
//            x.pel_id = NULL;
//            x.ppart = NULL;
//            x.pos_ = 0;
//            x.padjncy = NULL;
//            x.padjwgt = NULL;
        }

        CsrNodeInternal& operator=(const CsrNodeInternal &x){
            //*pxadj = *x.pxadj; // this is NOT changing
			*pvwgt = *x.pvwgt;
			*pel_id = *x.pel_id;
			*ppart = *x.ppart;
			pos_ = x.pos_;
			const int nx= *(x.pxadj+1) - *x.pxadj;
			assert(nx ==  (*(pxadj+1) - *pxadj) );
			//if not - hm... we have a problem....
            std::copy(x.padjncy,x.padjncy + nx, padjncy);
            std::copy(x.padjwgt,x.padjwgt + nx, padjwgt);
            std::copy(x.padj_neig_no,x.padj_neig_no + nx, padj_neig_no);
            return *this;
		}

        CsrNodeInternal& operator++() {
            ++pxadj;
            ++pvwgt;
            ++pel_id;
            padjwgt += *pxadj - *(pxadj-1);
            padjncy += *pxadj - *(pxadj-1);
            ++ppart;
            ++pos_;
            padj_neig_no += *pxadj - *(pxadj-1);
            return *this;
        }

        CsrNodeInternal& operator+(int n) {
            return move_forward(n);
        }

        CsrNodeInternal& operator--() {
            --pxadj;
            --pvwgt;
            --pel_id;
            padjwgt -= *(pxadj+1) - *pxadj;
            padjncy -= *(pxadj+1) - *pxadj;
            --ppart;
            --pos_;
            padj_neig_no -= *(pxadj+1) - *pxadj;
            return *this;
        }

        CsrNodeInternal& operator-(int n) {
            return move_back(n);
        }

        std::ptrdiff_t operator-(const CsrNodeInternal & x) const {
            return pel_id - x.pel_id;
        }

        bool    operator == (const CsrNodeInternal & x) const {
            return pel_id == x.pel_id;
        }

        bool    operator != (const CsrNodeInternal & x) const {
            return pel_id != x.pel_id;
        }

        bool operator < (const CsrNodeInternal & x) const {
            return pel_id < x.pel_id;
        }

        CsrNodeInternal& move_forward(const int n) {
            pxadj+=n;
            pvwgt+=n;
            pel_id+=n;
            padjwgt += *pxadj - *(pxadj-n);
            padjncy += *pxadj - *(pxadj-n);
            ppart+=n;
            pos_+=n;
            padj_neig_no += *pxadj - *(pxadj-n);
            return *this;
        }

        CsrNodeInternal& move_back(const int n) {
            pxadj-=n;
            pvwgt-=n;
            pel_id-=n;
            padjwgt -= *(pxadj+n) -*pxadj;
            padjncy -= *(pxadj+n) -*pxadj;
            ppart-=n;
            pos_-=n;
            padj_neig_no -= *(pxadj+n) -*pxadj;
            return *this;
        }

        Tind pos() const {return pos_;}
        Tind part() const {return *ppart; }
//    protected:
        pTind pxadj,padjncy,pvwgt,padjwgt,pel_id,ppart, padj_neig_no;
        Tind pos_;
    };


//////////////////////////////////////////////////////////////
class iterator
            : public std::iterator< std::random_access_iterator_tag, CsrNodeInternal  >
    {
    public:
        iterator(CSR *container = NULL)
            : node(container){}

        iterator(const iterator &x) {
            *this = x;
        }

        iterator( iterator &x) {
            node.pxadj = x.node.pxadj;
            node.pvwgt = x.node.pvwgt;
            node.pel_id = x.node.pel_id;
            node.ppart = x.node.ppart;
            node.pos_ = x.node.pos_;
            node.padjncy = x.node.padjncy;
            node.padjwgt = x.node.padjwgt;
            node.padj_neig_no = x .node.padj_neig_no;
            x.node.pxadj = NULL;
            x.node.pvwgt = NULL;
            x.node.pel_id = NULL;
            x.node.ppart = NULL;
            x.node.pos_ = 0;
            x.node.padjncy = NULL;
            x.node.padjwgt = NULL;
            x.node.padj_neig_no = NULL;
        }

        iterator& operator=(const iterator &it) {
            node.pxadj = it.node.pxadj;
            node.pvwgt = it.node.pvwgt;
            node.pel_id = it.node.pel_id;
            node.padjwgt = it.node.padjwgt;
            node.padjncy = it.node.padjncy;
            node.ppart = it.node.ppart;
            node.pos_ = it.node.pos_;
            node.padj_neig_no = it.node.padj_neig_no;

            return *this;
        }

        iterator& operator++() {
            ++node;
            return *this;
        }

        iterator& operator+(const int n) {
            node.move_forward(n);
            return *this;
        }

        iterator& operator--() {
            --node;
            return *this;
        }

        iterator& operator-(const int n) {
            node.move_back(n);
            return *this;
        }

        std::ptrdiff_t operator-(const iterator & x) const {
            return node - x.node;
        }

        bool    operator == (const iterator & x) const {
            return node == x.node;
        }

        bool    operator != (const iterator & x) const {
            return node != x.node;
        }

        bool operator < (const iterator & x) const {
            return node < x.node;
        }

        CsrNodeInternal& operator*() {
            return node;
        }

    protected:
        CsrNodeInternal node;
    };
///////////////////////////////////////////////////////////////////////
    CSR(const Tind size=0) {
        resize(size);
    }

    void resize (Tind size) {
        if(size > 0) {
            xadj_.resize(size+1,0);
            adjncy_.resize(2*(size+1)*n_max_neig,0);
            vwgt_.resize(size,0);
            adjwgt_.resize((size+1)*n_max_neig,0);
            el_loc_id_.resize(size,0);
            part_.resize(size,0);
            adj_neig_no_.resize(adjncy_.size(),0);
        }
        else {
            clear();
        }
    }

    void reserve (Tind size) {
        if(size > 0) {
            xadj_.reserve(size+1);
            adjncy_.reserve(2*(size+1)*n_max_neig);
            vwgt_.reserve(size);
            adjwgt_.reserve(2*(size+1)*n_max_neig);
            el_loc_id_.reserve(size);
            part_.reserve(size);
            adj_neig_no_.reserve(adjncy_.capacity());
        }
    }

    Tind capacity() const {
        return el_loc_id_.capacity();
    }

    CsrNodeInternal& at(const int pos) {
        assert(pos >= 0);
        assert(pos < size());
        return CsrNodeInternal(this,pos);
    }

    void clear() {
        xadj_.clear();
        adjncy_.clear();
        vwgt_.clear();
        adjwgt_.clear();
        el_loc_id_.clear();
        part_.clear();
        adj_neig_no_.clear();
    }

    iterator begin()  { return iterator(this); }
    iterator end()  {return iterator(this)+el_loc_id_.size(); }

    //const_iterator begin() const  { return const_iterator(this); }
    //const_iterator end() const  {return const_iterator(this,xadj_.size()); }

    Tind size() const { return el_loc_id_.size() ;}
    bool empty() const { return el_loc_id_.empty(); }

    pTind xadj() {return xadj_.data(); }
    pTind adjncy() {return adjncy_.data(); }
    pTind vwgt() {return vwgt_.data(); }
    pTind adjwgt() {return adjwgt_.data(); }
    pTind el_id() {return el_loc_id_.data(); }
    pTind part() {return part_.data(); }
    pTind adj_neig_no() { return adj_neig_no_.data(); }

    cpTind xadj() const {return xadj_.data(); }
    cpTind adjncy() const  {return adjncy_.data(); }
    cpTind vwgt() const  {return vwgt_.data(); }
    cpTind adjwgt() const  {return adjwgt_.data(); }
    cpTind el_id() const {return el_loc_id_.data(); }
    cpTind part() const {return part_.data(); }
    cpTind adj_neig_no() const {return adj_neig_no_.data(); }

    void push_back(CsrNode & val) {
        push_back(val.el_id, val.neighs, val.n_neighs, val.vwgt, 0, val.neig_no);
    }

    void push_back(CsrNodeInternal & val) {
        push_back(*val.pel_id, val.padjncy, *(val.pxadj+1) - *val.pxadj, *val.pvwgt, *val.ppart, val.padj_neig_no);
    }

    void push_back(Tind el_id, pTind neigs,Tind n_neigs,Tind vtxwgt,Tind part, pTind neig_no) {
        assert(el_id > 0);
        assert(neigs != NULL);
        assert(n_neigs > 0);
        assert(n_neigs <= n_max_neig);
        assert(vtxwgt>0);

        if(xadj_.empty()) {
            xadj_.push_back(0); // initial elem
        }

        el_loc_id_.push_back(el_id);
        vwgt_.push_back(vtxwgt);
        adjncy_.insert(adjncy_.end(),neigs,neigs+n_neigs);
        adjwgt_.insert(adjwgt_.end(),n_neigs,vtxwgt);
        xadj_.push_back( adjncy_.size() );
        part_.push_back(0);
        adj_neig_no_.insert(adj_neig_no_.end(), neig_no, neig_no+n_neigs);
    }

    template<class Archive>
    void serialize(Archive & archive)
    {
      //archive( xadj_, adjncy_, vwgt_, adjwgt_, el_loc_id_, part_ );
      archive( xadj_);
      archive( adjncy_);
      archive( vwgt_ );
      archive( adjwgt_ );
      archive( el_loc_id_ );
      archive( part_ );
      archive( adj_neig_no_ );
    }

    std::vector<T>& getOverlapElems(std::vector<T> & overlap, const int part_no)
    {
        overlap.clear();

        assert(part_no >= 0);

        for(iterator it = begin(); it != end(); ++it) {
            if(* (*it).ppart == part_no) { // if element belongs to part_no
                const int n_neigs = *((*it).pxadj+1) - *((*it).pxadj);
                for(int n=0 ; n <  n_neigs; ++n) {
                    const int neig_el = (*it).padjncy[n];
                    if(part_ [ neig_el ] != part_no) { // if has neighs from other part
                        overlap.push_back( el_loc_id_[neig_el] ); // add to overlap
                    }
                }
            }
        }

        std::sort(overlap.begin(), overlap.end() );
        overlap.erase( std::unique( overlap.begin(), overlap.end() ), overlap.end() );

        return overlap;
    }

//protected:
    std::vector<Tind> xadj_,adjncy_,vwgt_,adjwgt_,el_loc_id_,part_,adj_neig_no_;

};




//template<class Tind>
//bool operator!= (const typename CSR<Tind>::iterator & x, const typename CSR<Tind>::iterator& y) {
//    return x.node.pel_id != y.node.pel_id;
//}

//template<class Tind>
//bool operator<(const typename CSR<Tind>::iterator & x, const typename CSR<Tind>::iterator& y) {
//    return x.node.pel_id < y.node.pel_id;
//}

template<class Tind>
bool operator!= (const typename CSR<Tind>::CsrNodeInternal & x, const typename CSR<Tind>::CsrNodeInternal& y) {
    return x.pel_id != y.pel_id;
}

template<class Tind>
bool operator<(const typename CSR<Tind>::CsrNodeInternal & x, const typename CSR<Tind>::CsrNodeInternal & y) {
    return x.pel_id < y.pel_id;
}


#endif //_CSR_HPP_
