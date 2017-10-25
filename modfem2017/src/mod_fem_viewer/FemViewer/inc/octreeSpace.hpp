/*
 * octreeSpace.hpp
 *
 *  Created on: 03-09-2011
 *      Author: Paweł Macioł
 */

#ifndef OCTREESPACE_HPP_
#define OCTREESPACE_HPP_

#include <stack>
#include <algorithm>
#include "fv_float.h"
#include "MathHelper.h"
#include "BBox3D.h"

namespace FemViewer {

/* structure for nodes identyfication */

template<typename TObject,
	     typename TContent = ArrayT<TObject> >
class Octree
{
	typedef TObject 		 OctElement;
	typedef fvmath::CVec3<CoordType> OctPoint;
	typedef TContent 		 OctElems;
	typedef typename OctElems::iterator OctElemIter;
  protected:
	// Internal octree subnodes
	enum eNode { leaf, branch };
	// Forward declaration of internal iterator
	class iterator;
	// Base node
	class INode {
	public:
		virtual eNode type() const = 0;
		//virtual ~INode() {;}
	};

	// Leaf
	class Leaf : public INode {
	private:
		// Current level
		int			_resolution;
		// Extens of a leaf
		OctPoint 	_position;
		// List of elements in leaf
		OctElems	_elems;

	public:
		eNode type() const { return leaf; }
		      int& resolution()          { return _resolution; }
		const int& resolution()    const { return _resolution; }
		const OctPoint& position() const { return _position; }

		void addElem(const OctElement& el)
		{
			_elems.insert(el);
		}

		void deleteElem(const OctElement& el)
		{
			_elems.erase(el);
		}

		Leaf(const OctPoint& pos, int res = 1) : _resolution(res), _position(pos), _elems()
		{;}
		Leaf(const Leaf& lf) : _resolution(lf._resolution), _position(lf._position), _elems(lf._elems)
		{;}
		virtual ~Leaf()
		{;}
	};

	// Branch
	class Branch : public INode {
	private:
		// 0 - reserve for leaf 1 - 8 branches
		INode* _children[1 + 8];
	public:
		eNode type() const { return branch; }
		INode*& child(size_t idx)
		{
			return _children[idx];
		}

		iterator& serach(const OctPoint& X, const OctPoint& center, const CoordType size,
						 iterator& cell)
		{
			uint_t num(0);
			uint_t atb[3];
			for (uint_t i(0);i<3;++i) {
				atb[i] = (X[i] > center[i]);
				num += (atb[i] << (2-i));
			}

			cell.addHandle(this);
			cell.level() = num;

			if (_children[num]==0) {
				cell.level() = -1;
				cell++;
			} else {
				CoordType cor = 0.5*size;
				OctPoint newcenter(center);

				for (int i(0);i<3;++i) newcenter[i] += (atb[i]*2 -1)*cor;

				if ((*_children[num]).type() == branch) {
					(static_cast<Branch&>(*_children[num])).search(X,newcenter,cor,cell);
				} else {
					cell._leaf = static_cast<Leaf*>(_children[num]);
				}
			}
			return cell;
		}

		void addLeaf(const OctPoint& center,CoordType size,Leaf* lf)
		{
			int num(0);
			int atb[3];
			for (int i=0;i<3;++i) {
				atb[i] = ((*lf).position()[i] > center[i]);
				num += (atb[i] << (2-i));
			}

			// get the part of octree
			if (_children[num] == 0) {
				// If branch is empty add leaf
				_children[num] = lf;
			} else {
				CoordType cor = CoordType(0.5)*size;
				OctPoint newcenter(center);

				for (int i(0);i<3;++i) newcenter += (atb[i]*2 -1)*cor;

				if ((*_children[num]).type() == branch) {
					(static_cast<Branch&>(*_children[num])).addLeaf(newcenter,cor,lf);
				} else {
					Leaf* oldlf = static_cast<Leaf*>(_children[num]);
					Branch* br = new Branch();
					_children[num] = br;
					br->addLeaf(center,cor,lf);
					br->addLeaf(center,cor,oldlf);
				}
			}
		}

		void insert(const OctPoint& cent,
				    const CoordType& size,
				    const int res,
				    const int num,
				    const OctElems& el,
				    const BBox3D& itembb)
		{
			if (num==0) {
				if (_children[0]==0) {
					int lres = std::min(_maxres,res);
					_children[0] = new Leaf(cent,lres);
				}
				static_cast<Leaf*>(_children[0])->addElem(el);
			} else {
				if (_children[num]==0) {
					_children[num] = new Branch();
				}

				//br->insert
				int num = 0;
				int atb[8];
				for (int i=0;i<8;++i) {
					atb[i] = it.v[i] > center[i]);
					num += (atb[i] << (2-i));
				}
				CoordType hs = CoordType(0.5) * size;
				OctPoint nc(cent);
				for(int i=0;i<3;++i) nc[i] += (atb[i]*2 -1)*hs;

				Branch* br(0);
				if (_children[++num]==0) {
					br = new Branch();
					_children[num] = br;
				} else {
					br = static_cast<Branch*>(_children[num]);
				}
				int lres = getResolution(el.bbox.maxDim() * 0.5);
				br->insert(nc,hs,lres,num,el);

			}
		}


		OctElemIter begin()
		{
			OctElemIter it;
			it.addHandle(this);
			it++;
			return it;
		}
		explicit Branch(){ memset(_children,0x0,sizeof(_children)); }
		Branch(const Branch& br){;}
		virtual ~Branch() {;}
	private:
	};



public:
	Octree(const BBox3D& bbox,uint_t res)
	:  _pos(), _root(0), _maxres(0), _maxsize(0), _minsize(0)
	{
		assert(bbox.isInitialized());
		assert(res > 0);
		this->setResolution(res);
	}
	~Octree() {;}
public:
	class iterator {
	public:
		friend class Branch;
		typedef std::pair<int,Branch*> lHandle;
	private:
		Leaf* _leaf;
		std::stack<lHandle> _handles;
	public:
		size_t  getdepth() { return _handles.szie(); }
		Branch* getcurrfather()
		{
			if (_handles.size() >0) {
				return _handles.top().second;
			}
			else return 0;
		}
		iterator& operator++(int)
		{
			_leaf = 0;
			Branch* fa = getcurrfather();
			if (fa != 0) {
				_handles.top().first++;
				if (_handles.top().first < 8) {
					INode* node = fa->child(_handles.top().first);
					if (node != 0) {
						if (node->type() == leaf) {
							_leaf = static_cast<Leaf*>(node);
						} else {
							Branch* br = static_cast<Branch*>(node);
							this->addHandle(br);
							(*this)++;
						}
					} else {
						(*this)++;
					}
				} else {
					_handles.pop();
					(*this)++;
				}
			}
			return(*this);
		}

		int& level() { return _handles.top().first; }
		const int& level() const { return _handles.top().first; }

		void addHandle(Branch* br)
		{
			lHandle h;
			h.first  = -1;
			h.second = br;
			_handles.push(h);
		}

		operator Leaf*() { return _leaf; }

		iterator() : _leaf(0) {;}
		iterator(const iterator& it) : _leaf(it), _handles(it) {;}
		~iterator() {;}

	};

	iterator begin() const
	{
		iterator it;
		if (_root != 0) {
			if ((*_root).type() == leaf) {

			} else {
				Branch& br = static_cast<Branch&>(*_root);
				it = br.begin();
			}
		}
		return it;
	}

	iterator end() const
	{
		return iterator();
	}

	iterator search(const OctPoint& X)
	{
		iterator it;
		if (_root==0) {
			return this->end();
		} else if ((*_root).type() != branch) {
			it._leaf = static_cast<Leaf*>(_root);
			return it;
		} else {
			return (static_cast<Branch&>(*_root)).serach(X,_bbmesh.getCenter(),size,it);
		}

		return it;
	}

	void setResolution(const uint_t res)
	{
		if (!(res && (res & (res-1)))) _maxres = fvmath::npow2(res);
		else _maxres = res;
	}

	void setDimensions(const BBox3D& bbox)
	{
		if (_maxres > 0)
		{
			_pos = bbox.center();
			_maxsize = bbox.maxDim();
			_minsize = _maxsize / _maxres;
		}
	}

	void add(const OctElement& c, const BBox3D& elbb)
	{
		// If octree was not inited
		if (_root==0) {
			_root = new Leaf(elbbo.center());
			static_cast<Leaf&>(*_root).addElem(c);
		} else {
			if ((*_root).type() == leaf) {
				Leaf* lf = static_cast<Leaf*>(_root);
				/*if (_size==0) {
					setDimensions(lf->position(),position);
				}*/

				Branch* br = new Branch();
				_root = br;
				br->addLeaf(_pos,_size,lf);
				Leaf* nlf = new Leaf(position);
				nlf->addElem(c);
				br->addLeaf(_center,_size,nlf);
			} else {
				Leaf* nlf = new Leaf(position);
				nlf->addElem(c);
				(static_cast<Branch&>(*_root)).addLeaf(_center,_size,nlf);

			}
		}
	}

	void insert(const OctElement& item,const BBox3D& itembb)
	{
		OctPoint pos(_pos);
		int noct = check(pos,_maxsize*0.5,itembb);
		assert((0<=noct) && (noct <= 9));

		if (noct==0)
		{
			// Octree is empty
			if ( _root==0) {
				_root = new Leaf(_pos);
			}

			// There is only a leaf
			if ((*_root).type() == leaf) {
				(static_cast<Leaf*>(_root))->addElem(item);
			}
			else {
				Branch* br = static_cast<Branch*>(_root);
				Leaf* lf;
				if (br->child(0)==0)
				{
					lf = new Leaf(_pos,1);
					br->child(0) = lf;
				}
				else {
					lf = static_cast<Leaf*>(br->child(0));
				}
				lf->addElem(item);
			}

		}
		else {
			Branch* br;
			//Leaf* lf;
			if (_root==0) {
				//_root = new Branch();
				br = new Branch();
				static_cast<Branch*>(_root)->child(noct) = br;
			}

			if ((*_root).type() == leaf) {
				Leaf* lf = static_cast<Leaf*>(_root);
				br = new Branch();
				_root = br;
				br->child(0) = lf;
			} else {
				br = static_cast<Branch*>(_root);
			}

			br->insert(pos,_maxsize*0.5,_maxres,item,itembb);
		}
	}

	static int check(const OctPoint& pos, const CoordType size, const BBox3D& elem);
	inline BBox3D genBBoxForNode(const Vec3D& cent,int level)
	{
		assert(((level -1) & level) == 0);
		if (level==0) {
			return BBox3D(_pos,_maxsize*0.5f);
		} else {
			float size = this->_maxsize / static_cast<float>(1 << level);
			return BBox3D(cent,size);
		}
	}

	inline int getResolution(const T& size)
	{
		assert(size > 0);
		assert(_minsize > 0);
		int aln = npow2( static_cast<float>( ceil( (size / _minsize) ) );
		int lg2 = log2(aln);
		return this->_maxres >> lg2;
	}

private:
	OctPoint	 _pos;
	Branch*   	 _root;
	unsigned int _maxres;
	CoordType	 _maxsize;
	CoordType    _minsize;

};


template<typename TObject,typename TContent>
int
Octree<TObject,TContent>::check(const OctPoint& pos, const CoordType hs, const BBox3D& elbb)
{
	int num;
	int atb[8];
	OctPoint bc(elbb.center());

	atb[0] = (bc.x > pos.x);
	num = atb[0] << 2;
	atb[1] = (bc.y > pos.y);
	num += atb[1] << 1;
	atb[2] = (bc.z > pos.z);
	num += atb[2];

	//const CoordType lh = hs * 0.5;
	assert(hs > 0);
	OctPoint mv(hs,hs,hs);

	switch(num)
	{
	case 0:
		 mv.x *= -1;
		 mv.y *= -1;
		 mv.z *= -1;
		 break;
	case 1:
		mv.x *= -1;
		mv.y *= -1;
		break;
	case 2:
		mv.x *= -1;
		mv.z *= -1;
		break;
	case 3:
		mv.x *= -1
		break;
	case 4:
		mv.y *= -1;
		mv.z *= -1;
		break;
	case 5:
		mv.y *= -1;
		break;
	case 6:
		mv.z *= -1;
		break;
	case 7:
		break;
	}

	pos += mv;
	BBox3D lbb(pos,hs);
	unsigned int lo, up;
	lo = lbb.IsInside(elbb.mn) ? 1 : 0;
	up = lbb.IsInside(elbb.mx) ? 1 : 0;

	return (lo && up) ? num + 1 : 0;
}


} // end namespace FemViewer

#endif /* OCTREESPACE_HPP_ */


