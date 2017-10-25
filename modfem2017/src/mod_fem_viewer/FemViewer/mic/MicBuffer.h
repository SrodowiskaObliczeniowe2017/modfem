/*
 * MicBuffer.h
 *
 *  Created on: 20 lis 2015
 *      Author: dwg
 */

#ifndef MICBUFFER_H_
#define MICBUFFER_H_

namespace FemViewer {
namespace MIC {

class MICBuffer {
public:
	MICBuffer(){}
	virtual ~MICBuffer(){}

	void Map();
	void UnMap();

private:
	void*  m_data_ptr;
	size_t m_size;
};

} // end namespace MIC
} // end namespace FemViewer

void MICBuffer::Map()
{

}

void MICBuffer::UnMap()
{

}

#endif /* MICBUFFER_H_ */
