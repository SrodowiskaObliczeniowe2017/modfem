#include "fv_conv_utils.h"
#include <string.h>

bool IsDig(std::string &str)
{
    for(size_t i = 0; i < str.length(); ++i)
    {
        switch(str[i])
        {
            case '0': break;

            case '1': break;

            case '2': break;

            case '3': break;

            case '4': break;

            case '5': break;

            case '6': break;

            case '7': break;

            case '8': break;

            case '9': break;

            default : return false;
        }


    }

    return true;

}

bool IsCyfr(std::string &str)
{

    int x = 0;

    if( (str[0] == '+') || (str[0] == '-') )    
       ++x;

    //sprawdz ile jest kropek - je�li wiecej niz jedna to bl�d

    //je�li jest kropka to sprawdzic czy przed nia wystepuje liczba
    //lub sprawdzic czy:
    //nie jest na pierwszym miejscu
    //jesli na drugim to czy na pierszym jest liczba

    for(size_t i = x; i < str.length(); ++i)
    {

        //sprawdzenie czy znak jest od 0-9
        //jesli nie to return false

        switch(str[i])
        {
            case '0': break;

            case '1': break;

            case '2': break;

            case '3': break;

            case '4': break;

            case '5': break;

            case '6': break;

            case '7': break;

            case '8': break;

            case '9': break;

            case '.': break;

            default : return false;
        }

    }//for(unsigned int i = 0; i < str.length(); ++i)

    return true;


}

void int2bin(unsigned int nr,unsigned int size,char buff[])
{

   memset(buff,(int)'0',size);
   buff[--size] = '\0';

   do {
	   if (nr & 0x01) buff[size-1] = '1';
	   nr >>= 1;
   } while (--size);
}

void  int2bin2(unsigned int n, unsigned int  buffer_size, char  *buffer)
{
        unsigned int i = (buffer_size - 1);
        buffer[i] = '\0';

        while(i > 0)
        {
                if(n & 0x01)
                {
                        buffer[--i] = '1';
                } else
                {
                        buffer[--i] = '0';
               }

                n >>= 1;
        }
}
