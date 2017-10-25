#ifndef _LEGEND_H_
#define _LEGEND_H_

#include <string>
#include <vector>
#include <memory>
#include "fv_config.h"
#include "LegendValue.h"
#include "MathHelper.h"
#include "Color.h"
#include "ColorBar.h"
#include "GraphicElem.hpp"

namespace FemViewer {


/* Type of color map */
enum LegendType
{
	// Default colored style
    SOLID	 = 0,
    GRADIENT = 1
};

template<typename T>
struct LegendData {
	fvmath::Vec3<T> color;
	int   numBrreaks;
	T	  startValue;
	T	  deltaValue;
};


class Legend
{
  private:
	// Default number of bounds of legend
	static const int _ncont = 20;
	// Default legend type
	static const LegendType _dftype = GRADIENT;
	// Handle to GL representation
	std::shared_ptr<ColorBar> _hBar;
	// Minimal and maximal values stored for legend
    double 	_min, _max;
    // Begin and end of values in [_min, _max] range
    double 	_start, _end;
    // Current number of ranges
    unsigned int _nranges;
    // Legend type
    LegendType _type;
    // Whether the color black is the color outside [_start,_end] range
    bool _outsideBlack;
    // Status of a legend
    bool _ok;
    // list of color ranges
    std::vector<colormap_value> _colors;

  public:
    // Current legend color and value
    colormap_value color;

  public:
    // Ctrs and dtr
    explicit Legend();
    explicit Legend(const Legend & legend);
    explicit Legend(const std::string& path);
    ~Legend();
    // Operators
    Legend & operator=(const Legend & legend);
    bool operator ==(const Legend& rhs) const;
    void Reset();
    // Accessors
    double GetMin() const  { return _min; }
    double GetMax() const  { return _max; }
    double GetStart() const { return _start; }
    double GetEnd() const  { return _end; }
    unsigned int GetRangesNumber(void) const { return _nranges;}
    LegendType GetLegendType(void) const { return _type; }
    bool GetOutsideBlack(void) const {return this->_outsideBlack; };
    bool IsValid() const { return _ok; }
    std::vector<colormap_value>& GetColors() { return _colors; }
    const std::vector<colormap_value>& GetColors() const { return _colors; }
    template<typename T> int PackValuesIntoArray(T buffer[][4],const int max_size) const;

    // Setters
    void SetMin(double minValue) {
    	if (_min != minValue)
    		Set(minValue,_max,_start,_end,_nranges,_type,_outsideBlack);
    }
    void SetMax(double maxValue) {
    	if (_max != maxValue)
    		Set(_min,maxValue,_start,_end,_nranges,_type,_outsideBlack);
    }
    void SetStart(double startValue) {
    	if (_start != startValue)
       	Set(_min,_max,startValue,_end,_nranges,_type,_outsideBlack);
    }
    void SetEnd(double endValue) {
    	if (_end != endValue)
    		Set(_min,_max,_start,endValue,_nranges,_type,_outsideBlack);
    }
    void SetRange(unsigned int numRanges) {
    	if (_nranges != numRanges)
    		Set(_min,_max,_start,_end,numRanges,_type,_outsideBlack);
    }
    void SetType(LegendType type) {
    	if (_type != type)
    		Set(_min,_max,_start,_end,_nranges,type,_outsideBlack);
    }
    void SetOutsideBlack(bool outsideBlack) {
    	if (_outsideBlack != outsideBlack)
    		Set(_min,_max,_start,_end,_nranges,_type,outsideBlack);
    }

    void Draw(void)   { assert(_hBar); _hBar->Draw(); }
    void Update(void) { assert(_hBar); _hBar->Create(this); }
    void Create(void);
    // Set legend parameters
    void Set(double minValue = 0.0,
    		 double maxValue = 1.0,
    		 double startValue = 0.0,
			 double endValue = 1.0,
			 unsigned int rangesNumber = _ncont,
			 LegendType legendType = _dftype,
			 bool outsideBlack = false);


	//unsigned int GetBordersCount(void)
	//! Zwraca liczbe granic w legendzie nowa wersja
	int GetLiczbaGranic(void) const
	{
		if (_type == SOLID) {
			unsigned int a = _colors.size() - 1;
			a >>= 1; // /= 2;
			return ++a;
		}
		else if (_type == GRADIENT) {
			return (int)_colors.size();
		}

		return(-1);
	}

	//! Zwraca wartosc granicy opodanym indeksie 0..n-1
	double GetWartoscGranicy(const unsigned int &i) const
	{

		if (_type == SOLID) {
			return _colors[(i << 1)].value;
		}

		return _colors[i].value;
	}



        //! Zwraca numer granicy dla podnaje wartosci lub -1 jesli nie jest granica
        int IsBorderNr(double value) const;



        //! Pobiera wartosci dla kt�rych zwraca kolor�w je�li mieszcz� si� w tym samym przedziale legendy
        std::vector<ColorRGB> GetColors(const std::vector<double> & vValues);

        //! Pobiera element graficzny i wylicza dla niego wartosci kolorow legendy
        /*!
            \return zwraca false w przypadku gdy nie mozna wszytskich wartosci elementu zakwalifikowac
            do jednego przedzialu przez co nie mozna dla wsyztskich ustawic koloru
        */
        //bool GetColors(GraphElement2<double> & ge2) const;

        //podajemy wektor graficznyc elementow i dla kazego elementuwektora wykonywana jest
        //opcja podawania koloru
        //! Pobiera wektor element�w graficznych i wylicza dla niego wartosci kolor�w legendy
        //bool GetColors(std::vector< GraphElement2<double> > & ge2) const;

        //! Dla podanej wartosci zwraca mozliwe przedzialy w jakich sie znajduje
        /*!
            Dla granicy powinny to byc dwie warto�ci: lewa i prawa
            dla wartosci miedzy granicami jest to tylko jeden przedzial
        */
        //void GetRanges(const double & d, std::vector<unsigned int> & vUi);
        std::vector<unsigned int> GetRanges(const double & d) const;

        //! Dla podanych dwu wartosci podaje mozliwy przedzial lub -1 jesli niemoze ich zakwalifikowac do jednego luba zadnego
        int GetRanges(const double & v1, const double & v2) const;

        //! Dla przedzia�u o podanym numerze pobiera jego graniczne wartosci
        void GetRangeValues(unsigned int i, double & v1, double & v2) const;

        //! Zwraca liczbe przedzialow w legendzie
        size_t GetRangesSize(void) const;


        //! Dla podanego przedzialu i wartosci zwrocenie skladowe koloru
        /*!
            \return false - jesli wartosc jest z po za podanergo przedzialu,
            lub z jakiegos innego powodu nie mozna ustalic skladowych koloru dla podanych parametrow
        */
        bool GetColor(const double & value, const unsigned int & range, ColorRGB & c) const;

    /*
        //zwraca 0 dla zera, 1 dla koloru
        unsigned int GetType(void)
        {

            switch(this->legendType)
            {
                case COLOR:
                    return 0;
                case GRADIENT:
                    return 1;
            }

        }
    */


    public:
        void GetRange3(const unsigned int & range, double & sValue, ColorRGB & sColor, double & eValue, ColorRGB & eColor);
    protected:
        void GetRangeColor(const unsigned int & range, double & sValue, ColorRGB & sColor, double & eValue, ColorRGB & eColor);
        void GetRangeGradient(const unsigned int & range, double & sValue, ColorRGB & sColor, double & eValue, ColorRGB & eColor);
    public:

        //! Dla podanego przedzia�u zwraca poczatkowa oraz koncowa - wartosc oraz kolor
        /*!

            Je�li podany przedzialu jest z po za zakrestu to w wartosciach zwracanych nie jest dokonywana zmiana

            \param range numer sprawdzanego przedzialu
            \param sValue poczatkowa wartosc
            \param sColor poczatkowy kolor
            \param eValue koncowa wartos
            \param eColor koncowy kolor
        */
        void GetRange(const unsigned int & range, double & sValue, ColorRGB & sColor, double & eValue, ColorRGB & eColor)
        {
            if (range >= this->_colors.size())
                return;


            unsigned int rangeStart, rangeEnd;

            rangeStart = range;
            rangeEnd = range+1;

            sValue = _colors[rangeStart].value;
            eValue = _colors[rangeEnd].value;

            if ( (_colors[rangeStart].type & INSIDE) && (_colors[rangeEnd].type & INSIDE)   )
            {

                #ifdef DEBUG2
                    cout << endl << "*IN*" << endl;
                #endif
                //jesli obie granice s a wewnetrzne


                //w zaleznosci od ustawienia rodzaju legendy
                if (_type == SOLID)
                {
                    //jednolitym - czyli wartosc koloru pierwszego ??
                    sColor = _colors[rangeStart].GetColor();
                    eColor = _colors[rangeStart].GetColor();


                }
                else if  (_type == GRADIENT)
                {
                    //gradientem-przejsciem od do
                    sColor = _colors[rangeStart].GetColor();
                    eColor = _colors[rangeEnd].GetColor();


                }


            }
            else if ( (_colors[rangeStart].type & OUTSIDE) && (_colors[rangeEnd].type & OUTSIDE)   )
            {
                #ifdef DEBUG2
                    cout << endl << "*OUT*" << endl;
                #endif
                //jesli obie granice sa zewnetrzne

                //jesli wlaczony czarny to
                if (_outsideBlack)
                {
                    //jesli czarny
                    sColor = ColorRGB();
                    eColor = ColorRGB();
                }
                else
                {
                    //jesli wlaczony kolor to

                    //nalezy sprawdzic czy obie wartosci sa mniejszcze od start czy wieksze od end
                        //jesli mniejsze to kolor start
                        //jesli wieksze to kolor end
                    if ( ( _colors[rangeStart].value < _start ) && ( _colors[rangeEnd].value < _start) )
                    {
                        int i = GetStartIndex();
                        sColor = _colors[i].GetColor();
                        eColor = _colors[i].GetColor();

                    }
                    else if ( ( _colors[rangeStart].value > _end ) && ( _colors[rangeEnd].value > _end) )
                    {
                        int i = GetEndIndex();
                        sColor = _colors[i].GetColor();
                        eColor = _colors[i].GetColor();
                    }
                }
            }
            else
            {
                #ifdef DEBUG2
                    cout << endl << "*INOUT*" << endl;
                #endif
                //ostatni przypadek jesli jedna wew a druga zew
                if (_outsideBlack)
                {
                    //jesli czarny
                    sColor = ColorRGB(0.0, 0.0, 0.0);
                    eColor = ColorRGB(0.0, 0.0, 0.0);
                }
                else
                {
                    //jesli wlaczony kolor to

                    //sprawdzeie ktory jest inside i pobranie jego kolrow
                    if (_colors[rangeStart].type & INSIDE)
                    {
                        sColor = _colors[rangeStart].GetColor();
                        eColor = _colors[rangeStart].GetColor();
                    }
                    else if (_colors[rangeEnd].type & INSIDE)
                    {
                        sColor = _colors[rangeEnd].GetColor();
                        eColor = _colors[rangeEnd].GetColor();
                    }

                }//if (outsideBlack)
            }//else
        };

        //z sugerownaniem przedzialu,
        //jesli jest granica i na do wyboru dwie wartosci
        //a jedna z nich jest zgodna z sugerowana to wtedy wybiera sugerowana
        //jesli przedzial gugerowany jest niezgodny z zadna z wartosci to wybiera pierwsza


        /*! \fn double GetInterp(double value)
        *   \brief A member function.
        *
        *   Dla podanej wartosci rozwiazania zwraca wartosc interpolacji <0-1>
        *   z przynaleznego przedzialu tej wartosci
        *
        *   \param value sprawdzana wartosc.
        *   \return interpolacja dla podanej wartosci z przynale�nego zakresu
        */
        double GetInterp(double value);

        /*! \fn double GetInterp(double a, double b, double value)
        *   \brief A member function.
        *
        *   Dla podanej wartosci rozwiazania zwraca wartosc interpolacji <0-1>
        *   z podanego przedzialu od a do b
        *
        *   \param a poczatek zakresu.
        *   \param b koiec zakresu.
        *   \param value sprawdzana wartosc.
        *   \return interpolacja dla podanej wartosci z daneo zakresu
        */
        double GetInterp(double a, double b, double value) const;

        //! Na podstawie poczatkoej i koncowej wartosci oraz interpolacji zwraca wartosc wewnetrzna
        double GetValue(double a, double b, double interp) const;






        //! Dla dwoch wartosci zwraca wektor ktory zawiera je jako poczatkowa i koncowa oraz wewnatrz wszytskie granice miedzy tymiwartosciami
        void GetBorders(const double &s, const double &e, std::vector<double> &vD) const
        {

            double dTmp;

            vD.clear();

            vD.push_back(s);

            //sprawdzenie czy sa sobie rowne
            if ( fvmath::Compare(s,e) )
            {
                vD.push_back(e);
                return;
            }

            //sprawdzenie czy s i e sa w tej samej granicy
            if ( GetRanges(s, e) >= 0  )
            {
                vD.push_back(e);
                return;
            }

            if (s < e)
            {
                //#ifdef DEBUG4
                //    cout << "<" << endl;
                //#endif
                //vD.push_back(s);
                for (int i = 0; i < GetLiczbaGranic(); ++i)
                {
                    dTmp = GetWartoscGranicy(i);

                    if (  (dTmp > s) && (dTmp < e)  )
                        vD.push_back(dTmp);
                }
                //vD.push_back(e);
            }
            else if (s > e)
            {
                //#ifdef DEBUG4
                //    cout << ">" << endl;
                //#endif
                //vD.push_back(s);

                //for (unsigned int i = 0; i < GetLiczbaGranic(); ++i)
                for (int i = GetLiczbaGranic(); i > 0; --i)
                {
                    dTmp = GetWartoscGranicy(i-1);

                    //if (  (dTmp > e) && (dTmp < s)  )
                    //if (  (dTmp < s) && (dTmp > e)  )

                    //#ifdef DEBUG4
                    //    cout << dTmp << "<" << s << " && " << dTmp<< ">"<< e ;
                    //#endif

                    //if (  (s > dTmp) && (dTmp > e)  )
                    if (  (dTmp < s) && (dTmp > e)  )
                    {
                        vD.push_back(dTmp);
                        //#ifdef DEBUG4
                        //    cout << " tak";
                        //#endif
                    }
                    //cout << endl;
                }

                //vD.push_back(e);

            }
            vD.push_back(e);

        };









	protected:



    public:
        //! Dla podanej wartosci zwraca gorna oraz dolna wartosc przedzialu - nowa wersja
        bool GetRangesBorder(const double & value, double & s, double &e) const;

        //! sprawdza czy podana wartosc jest granica - nowa wersja
        bool IsBorder2(double value) const;

        bool IsBorder3(double value) const;

	protected:

        //! Zwraca pozycje min i max wzgledem start end
        /*!
            mi dla min, mx dla max:
            -1 - mnijeszy od min
            0 - miedzy min a max
            1 wiekszy od max
        */
        void MinMaxPos(int & mi, int & mx);

        //! Na podstawie wewnetrznych ustawien
        void FillStartEndColor(void);

        void FillStartEndGradient(void);

        void AddMinMaxColor(void);

        void AddMinMaxGradient(void);


    //! \name Stare metdy juz najprawdopodobniej niewykorzystywane


    //protected:
    public:


            //sprawdzenie czy min i max sa w odpowiedzeniej kolejnisci i to samo dla start end

		//void SetMinMax2(double min, double max);

        //ustawienia domysle z parametrami min, max
        //jaki start, end ustawia odpowiedni min, max
        //z podzialem na ilosc elementow i ustaleniem na tryb kolor lub gradient

        //! Dla podanych wartosci tworzy wewn�trzne strukt�re legendy.Ustawia wszytskie dane niezbedne do funkcjonowanie legendy.
        /*!
            Ustawia wszytskie dane niezbedne do funkcjonowanie legendy - powinna byc stosowana jako pierwsza funkcja
            po stworzeniu a przed faktycznym wykorzystaniem poniewaz wszytskie metody modyfikujace opieraja sie na niej
            zminiajac tylko niezb�dne paramtery. w przypadku nie uzycia jej lub kolejno metod modyfikujacych wszystkie parametry
            legenda mo�e dzia�a� nipoprawnie

            \param min najmijesza warto�� rozwi�zania.
            \param max najwi�ksza warto�c rozwi�zania.
            \param start pocz�tkowa warto�� rozwi�zania.
            \param end ko�cowa warto�c rozwi�zania.
            \param rangesNumber liczba przedzia��w.
            \param legendType sposob prezentacji kolor�w w przedzia�ach.
            \param outsideBlack rodzaj koloru po za wybranym roziwazaniem a wewnatrz dopuszczalnych wartosci.

        */
        //void SetDefault(double min, double max, double start, double end, unsigned int rangesNumber, LegendType legendType, bool outsideBlack);


        //!Dla podanych watosci min i max dokonuje podzialu na 5 progow (4 przedzialy). Dla ustalonych warto�ci kolor�w.
		//void SetMinMax
		//(
		//!min powiniem byc wiekszy od max
		//double min,
		//!max powinien byc wiekszy od min
		//double max
		//);

                /*! \fn int SetMinMax(double min, double max, double start, double end, unsigned int num)
        *   \brief A member function.
        *
        *   Sprawdzenie czy podana wartosc jest granica.
        *   Zwraca wartosc < 0 je�li blad parametrow wejsciowych lub instnieje w nich jakis konflikt
        *   W przypadku bledu wektor wartosci legendy pozostaje niezmieniony
        *   -1  - nieznany blad
        *   -1  - jesli num jest zla wartoscia
        *       - jesli poczatek i koniec jest z po za wartosci min i max
        *   ... - ...
        *
        *   \param min
        *   \param max
        *   \param start
        *   \param end
        *   \param num
        *   \return zwraca informacje na temat poprawnosci parametro wejsciowych
        */
        //int SetMinMax(double min, double max, double start, double end, unsigned int num);

        /*! \fn unsigned int GetIndex(double value)
        *   \brief A member function.
        *
        *   Zwraca indeks przedzialu <1-n> do ktorego nalezy na podstawie wartosci.
        *
        *   \param value sprawdzana wartosc.
        *   \return zwraca numer zakresu
        */
        //unsigned int GetIndex(double value);

        /*! \fn unsigned int GetIndexMax(void)
        *   \brief A member function.
        *
        *   //Zwraca indeks przedzialu <1-n> do ktorego nalezy na podstawie wartosci.
        *   //zwraca maksymalna wartosc indeksu
        *   Zwraca liczba indeksow (0-n)
        *
        *   \return zwrca liczbe bedaca liczba indeksow
        */
        int GetIndexMax(void){return (int)_colors.size();};

		/*! \fn void GetColor(double value)
        *   \brief A member function.
        *
        *   zwraca wartosc skladowych kolorw dla podanej wartosci
		*   w postaci publicznego skladnika klasy
        *
        *   \param value wartosc dla korej podane sa skaldowe kolorow.
        */
		void GetColor(double value);

		// Get 3 byte color
		float GetColor3b(double value) {
			GetColor(value, color);
			return color.rgb.r + color.rgb.g*256.0 + color.rgb.b*256.0*256.0;
		}

		double GetValue(const ColorRGB& rgb);

        /*! \fn void GetColor(double value, LegendValue &lv)
        *   \brief A member function.
        *
        *   Dla podanej wartosci zwraca w postaci drugiego parametru
        *   informacje na temat wartosci oraz skladowych kolorow jej odpowiadjacych
        *
        *   \param value sprawdzana wartosc.
        *   \param lv ????
        */
		void GetColor(double value, colormap_value &lv);

		        /*! \fn bool IsBorder(double value)
        *   \brief A member function.
        *
        *   Sprawdzenie czy podana wartosc jest granica - stara wersja.
        *
        *   \param value sprawdzana wartosc.
        *   \return prawda jesli wartosc jest granica, w przeciwnym razie false
        */
        bool IsBorder(double value);



        /*! \fn bool IsInStartEnd(double value)
        *   \brief A member function.
        *
        *   Sprawdza czy podana wartosc jest >= start i <=end. Zwraca true dla prawdy lub false jesli nie
        *
        *   \param value sprawdzana wartosc.
        *   \return prawda jesli wartosc jest granica, w przeciwnym razie false
        */
        bool IsInStartEnd(double value);


                //! Zwraca indeks elementu odpowiadajacego typu start w legendzie
        int GetStartIndex(void)
        {
            for (unsigned int i = 0; i < _colors.size(); ++i)
            {
                if (_colors[i].type ==  START)
                    return i;
            }
            return -1;
        }

        //! Zwraca indeks elementu odpowiadajacego typu end w legendzie
        int GetEndIndex(void)
        {

            for (unsigned int i = 0; i < _colors.size(); ++i)
            {
                if (_colors[i].type ==  END)
                    return i;
            }
            return -1;

        };

    public:

                /*! \fn unsigned int GetIndex(double value)
        *   \brief A member function.
        *
        *   1) Dla podanego zakresu koloru w postaci dwu liczb
		*   zwraca w postaci: min, progiwystepujace miedzy min a max, max.
		*   Jako skladowe wartosc i odpowiednio przyapsowane kolory.
		*   2) Dla podanej wartosci poczatkowej i koncowej zwraca wynik w posaci publiczbego wektora colorRange.
		*   Sa to kolejno sk�adowe (wartosc+skladowe koloro): poczatkowa, przedzialy jesli takowe wsyepujamiedzy wartosciami, koncowa
        *
        *   \param valueStart poczatkowa wartosc z zakresu min-max
        *   \param valueEnd konczowa wartosc z zakresu min-max
        */
        void GetInterval(double valueStart, double valueEnd);

#ifdef DEBUG
		void Dump()
		{
			std::cout	<< "Legend"
						<< "\n\tColor:\t\t\t" << color.AsString()
						<< "\n\tType:\t\t\t" << (_type ? "GRADIENT" : "SOLID")
						<< "\n\tMin/Max value:\t\t" << _min << " | " << _max
						<< "\n\tStart/End value:\t" << _start << " | " << _end
						<< "\n\t# of ranges\t\t" << _nranges
						<< "\n\tIs ou black:\t\t" << (_outsideBlack ? "TRUE" : "FALSE")
						<< std::endl;


		}
#endif
    public:


};

template<typename T>
int Legend::PackValuesIntoArray(T buffer[][4],const int max_size) const
{
	if (buffer == nullptr) return(-1);

	int num = static_cast<int>(_colors.size());
	num = num > max_size ? max_size : num;
	for (int i = 0; i < num; ++i) {
		buffer[i][0] = _colors[i].rgb.r;
		buffer[i][1] = _colors[i].rgb.g;
		buffer[i][2] = _colors[i].rgb.b;
		buffer[i][3] = static_cast<float>(_colors[i].value);
		//buffer[i] = static_cast<float>(_colors[i].value);
		//printf("value%d = %f\n",i,buffer[i][3]);
	}
	return num;
}




} // end namespace FemViewer

#endif /* _LEGEND_H_
  */
