#ifndef _FV_MATH_PARSER_H_
#define _FV_MATH_PARSER_H_

/* dictionary */
//#include"fv_dictstr.h"
#include<vector>
#include<stack>
#include<string>
#include<cstdlib>


namespace fvmathParser {
extern double my_apr_elem_calc_3D_mod(
		/* returns: Jacobian determinant at a point, either for */
		/* 	volume integration if Vec_norm==NULL,  */
		/* 	or for surface integration otherwise */
		int Control,	    /* in: control parameter (what to compute): */
				    /*	1  - shape functions and values */
				    /*	2  - derivatives and jacobian */
				    /* 	>2 - computations on the (Control-2)-th */
				    /*	     element's face */
		int Nreq,	    /* in: number of equations */
		int *Pdeg_vec,	    /* in: element degree of polynomial */
		int Base_type,	    /* in: type of basis functions: */
				    /* 	(APC_BASE_TENSOR_DG) - tensor product */
				    /* 	(APC_BASE_COMPLETE_DG) - complete polynomials */
	        double *Eta,	    /* in: local coordinates of the input point */
		double *Node_coor,  /* in: array of coordinates of vertices of element */
		double *Sol_dofs,   /* in: array of element' dofs */
		double *Base_phi,   /* out: basis functions */
		double *Base_dphix, /* out: x-derivatives of basis functions */
		double *Base_dphiy, /* out: y-derivatives of basis functions */
		double *Base_dphiz, /* out: z-derivatives of basis functions */
		double *Xcoor,	    /* out: global coordinates of the point*/
		double *Sol,        /* out: solution at the point */
		double *Dsolx,      /* out: derivatives of solution at the point */
		double *Dsoly,      /* out: derivatives of solution at the point */
		double *Dsolz,      /* out: derivatives of solution at the point */
		double *Vec_nor     /* out: outward unit vector normal to the face */
		);

//namespace fvmathParser {

	// Types of math elements
	enum ElemType {
		Number,
		Vector,
		Plus,
		Minus,
		Multiply,
		Divide,
		Power,
		Open,
		Close
	};

	// MathElement definition
	struct MathElement { 
		ElemType type; double value; 
	};


	class MathCalculator {
	public:

		//zwraca gtowy element na podstawie tekstu
		//1     - rozpoznano poprawnie
		//<0    - blad
		static int GetElement(ElemType &type, double &value, std::string &str);

        //zwraca wypelniony element
        //zapis elementu
        //maksymalny rozmiar wektora ktory nie moze byc pzekroczony podczas analizy, jesli wiekszy to blad
		static void GetElement(MathElement & el, std::string &str);



        //ciag funkcji
        //rozmiar wektora n (0:n-1)

        //static int Test(std::string func, int vectorSize);

		//testowanie poprawnosci sk��dni
		//analizuje

		//analizuje
		//zapisuje do vectora
		//wykonuje wyliczenia przechodzac po wektorze w kolejnosci wykoanywania dzia�a�
		//nie zmienia� pierwotnego wektora i generowanie nast�pnych na bazie poprzednich

        //-1 - nieznany lb koncowy blad
        //-2 - blad analizy elementu z funkcji wejsiowej
        //-3 - blad zakresu podania wektora
        //<-10  - blad elementu analizy n-tego elementu (liczac od zera) -10-n
		//static int ONP(std::string func, int vectorSize);
		static int ONP(std::string func, unsigned int vectorSize, std::vector<MathElement> &wyjscie);

        //podaje priorytet w zalezno�ci od rodzaju elementu
        //3     - ^     - POW
        //2     - *, /  - Multiply, Divide
        //1     - +, -  - PLUS, MIN
        //0     - (, )  - OPEN, CLOSE
        //-1    -       - NUM, VEC
        static int GetPriority(ElemType ElType);


        //zwraca wyjscie onp i wierzcho�ek stosu jako string
        //static std::string GetString(std::vector<Element> &wyjscie);
        //static std::string GetString(std::vector<Element> &vElement);
        static std::string GetString(std::vector<MathElement> &vElement, std::stack<MathElement> &stos);

        //zwraca wektor (wyjscie) onp jako string
        static std::string GetString(std::vector<MathElement> &vElement);

        //(2+3)*5
        //2 3 + 5 *

        //((2+7)/3+(14-3)*4)/2
        //2 7 + 3 / 14 3 - 4 * + 2 /

        //obliczenie wyrazenia onp
        //wykorzystanie Calculate
        // 1 - wszytsko wporzadku
        // <0 - blad
        static int ONPCalculate(std::vector<MathElement> &vElement, std::vector<double> v, double &result);

        //matematyczna operacja typu v
        //podawane sa dwa argumenty i znak opracji
        //jednym z argumentow lub oba moga by� wektorem
        //trzy elementy, wektor rozwiazania
        //dwa pierwsze sa wartoscia lub reprezentacja
        //trzeci znakiem
        //czwarty - wynik rozwiazania double
        //wektor rozwiazan - wektor typu double
        static void Calculate(MathElement &a, MathElement &b, MathElement &sign, MathElement &result, std::vector<double> v);

	};
} // end Namespace

#endif /* _FV_MATH_PARSER_H_
*/
