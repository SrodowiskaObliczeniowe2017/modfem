#ifndef _FV_STRUTLS_H_
#define _FV_STRUTLS_H_

#include<string>
#include<vector>

std::string trimString(const std::string& sText,
					   const std::string& sToken);

int trimString(const std::string& sText,
			   std::vector<std::string>& vTokens,
			   const std::string& sSeps = " \t");

int countWords(const std::string& sText);

int countSubString(const std::string& sText, 
				  const std::string& sToken);

int getSubStrPosition(const std::string& sText,
					  const std::string& stoken,
					  int num = 1);
#endif /* _FV_STRUTLS_H_
*/
