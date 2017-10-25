#include "fv_str_utils.h"
#include "fv_compiler.h"


std::string trimString(const std::string& sText,
					   const std::string& sToken)
{
	 std::string sNewText;

	  if (!sText.empty()) {
		std::string::size_type beg = sText.find_first_not_of(sToken);
		std::string::size_type end = sText.find_last_not_of(sToken);

		if (beg != std::string::npos && end != std::string::npos) {
		  sNewText = sText.substr(beg, (end+1) - beg);
		}
		else {
		  assert((beg == std::string::npos) && (end == std::string::npos));
		}
	  }

	  return sNewText;

}

int trimString(const std::string& sText, std::vector<std::string>& vTokens,
			   const std::string& sSeps)
{
	// Clear vector
	vTokens.clear();
	// Skip delimiters at beginning.
	std::string::size_type lpos = sText.find_first_not_of(sSeps, 0);
	// Find first "non-delimiter".
	std::string::size_type pos     = sText.find_first_of(sSeps, lpos);

	while (std::string::npos != pos || std::string::npos != lpos)
	{
		// Found a token, add it to the vector.
		vTokens.push_back(sText.substr(lpos, pos - lpos));
		// Skip delimiters.  Note the "not_of"
		lpos = sText.find_first_not_of(sSeps, pos);
		// Find next "non-delimiter"
		pos = sText.find_first_of(sSeps, lpos);
	}

	return static_cast<int>(vTokens.size());
}


int countWords(const std::string& sText)
{
	int  cnt = 0;
	bool flag = false;
	std::string::size_type l = sText.size();
	
	for(std::string::size_type i=0; i<l; ++i) {

		if(sText[i] == ' ') {
		  flag = false;
		}
		else {
		  if(!flag) {
			++cnt;
		  }
		  flag = true;
		}
	  }
	  return cnt;

}

int countSubString(const std::string& sText,
				  const std::string& sToken)
{
	if(sText.empty() || sToken.empty())
		return(-1);

    int cnt = 0;
	std::string::size_type poz = 0, w;

    do {
		w = sText.find(sToken, poz);
		if(w != std::string::npos)
        {
			++cnt;
			poz = w + 1;
        }

	} while(w != std::string::npos);

    return( cnt);

}


int getSubStrPosition(const std::string& sText,
					  const std::string& sToken,
					  int num)
{
	if(sText.empty() || sToken.empty() || num < 1)
		return(-1);

	std::string::size_type poz = 0;

    for(int i = 0; i < num; ++i)
    {
        poz = sText.find(sToken, poz);
		if(poz == std::string::npos)
            return(-1);
        else
			poz++;
    }

    return( (int)--poz);
}
