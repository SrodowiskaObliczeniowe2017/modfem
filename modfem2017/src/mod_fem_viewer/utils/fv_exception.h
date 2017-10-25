#ifndef _FV_EXCEPTION_H_
#define _FV_EXCEPTION_H_

#include<stdexcept>
#include<string>
#include<ostream>


//namespace FemViewer
//{
	class fv_exception : std::exception
	{

	public:

		fv_exception(): sReason("Unknown exception") {}
		fv_exception(const char* resson) : sReason(resson) {}

		fv_exception(const std::string& reason)
			: sReason(reason) 
		{}

		fv_exception(const fv_exception& excpt) 
			: sReason(excpt.sReason) 
		{}

		fv_exception& operator =(const fv_exception& excpt)
		{
			if(this != &excpt)
				sReason = excpt.sReason;
			return *this;
		}

		virtual ~fv_exception() throw() {}

		virtual const char* what() const throw() { return sReason.c_str(); }

	protected:

		std::string sReason;

		friend std::ostream& operator <<(std::ostream &lhs, const fv_exception& excpt);
	};

	inline std::ostream& operator <<(std::ostream &lhs, const fv_exception &rhs)
	{
		return lhs << rhs.sReason;
	}
//}

#endif /* _FV_EXCEPTION_H_
*/
