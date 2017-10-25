#ifndef Handle_H_
#define Handle_H_


template< typename T >
class Handle  
{
public:
	Handle() : _handle(0) {}
	virtual ~Handle() { _handle = T(0); }
	operator T() { return _handle; }
	virtual int  IsInit(void) { return 0; }
	virtual int  Init(int argc,char** argv) { return 0; }
	virtual void ShutDown(void) {} 
protected:
	T _handle;
};

#endif /* Handle_H_
*/