#ifndef _PCSR_HPP_
#define _PCSR_HPP_

#include "Csr.hpp"

template <typename Tind>
class PCSR : public CSR<Tind>
{
public:

    friend std::ostream& operator<<(std::ostream & os,  PCSR<Tind> const& pcsr) {
        return pcsr.write(os);
    }

    friend std::istream& operator>>(std::istream & is, PCSR<Tind> & pcsr) {
        return pcsr.read(is);
    }

    typename std::vector<Tind>::iterator write(std::vector<Tind> & v) const {
        CSR<Tind>::write(v);

        v.push_back(vtxDist_.size());
        v.insert(v.end(),vtxDist_.begin(),vtxDist_.end());

        return v.end();
    }

    std::ostream& write(std::ostream & os) const {
        CSR<Tind>::write(os);

        os << vtxDist_.size() << " ";
        //for(const auto & el: vtxDist_) { os << el << " " ; }
        for(int i=0; i < vtxDist_.size(); ++i) { os << vtxDist_[i] << " "; }
        os << "\n\n";
        //os.write(reinterpret_cast<const char*>(vtxDist()),vtxDist_.size()*sizeof(Tind));
        return os;
    }

    typename std::vector<Tind>::const_iterator  read(const std::vector<Tind> & v) {
        typename std::vector<Tind>::const_iterator it = CSR<Tind>::read(v);

        vtxDist_.resize(*it); ++it;
        std::copy(it,it+vtxDist_.size(),vtxDist_.begin());

        return it;
    }

    std::istream& read(std::istream & is) {
        CSR<Tind>::read(is);

        Tind tmp_size = 0;
        is >> tmp_size;
        vtxDist_.resize(tmp_size);
        std::istream_iterator<Tind> iit(is);
        std::copy(iit, iit+vtxDist_.size(),vtxDist());
        //is.read(reinterpret_cast<char*>(vtxDist()),vtxDist_.size()*sizeof(Tind));
        return is;
    }

    PCSR(const int nProc = 0, const Tind size =0): CSR<Tind>(size) , vtxDist_(nProc,0)
    {}

    void resize_n_proc(const int size) {
        // vtxDist is +1 size (last number is total number of graph nodes)
        vtxDist_.resize(size+1,0);
    }

    int n_proc() const { return vtxDist_.size()-1 ;}

    Tind* vtxDist() { return vtxDist_.data(); }
    const Tind* vtxDist() const { return vtxDist_.data(); }

    void clear() {
        vtxDist_.clear();
        CSR<Tind>::clear();
    }

    bool    empty() const {
        return vtxDist_.empty();
    }

    template<class Archive>
    void serialize(Archive & archive)
    {
        CSR<Tind>::serialize(archive);
        archive( vtxDist_ );
    }

protected:
    std::vector<Tind> vtxDist_;
};



#endif //_CSR_HPP_
