FemViewr - aplikacja graficzna do wizualizacji rozwi¹zañ MES na siatkach pryzmatycznych i czworoœciennych
Autor: Pawe³ Macio³, Politechnika Krakowska
wxWidgtes-version: 2.9.1.
Kompilator: Microsoft VC++ 2008.
System testowy: Windows 7.

FemViewer dzia³a w oparciu o biblioteki wxWidgets, w tym celu nalezy mieæ skompilowane biblioteki 
wxWidgets w postaci dynamicznej (.dll w Windows).

Przygotowanie (Win 7):
1. pobraæ i zainstalowaæ najnowsz¹ wersje bibliotek wxWidgets dla Win7.
2. w zmienych systemowych ustawic zmien¹ WXWIN na scie¿kê\do\katalogu\wxWidgets (np. C:\wxWidgets-2.9.1)
3. wejœæ do katalogu sciezka\do\wxWidgets\include\wx\msw\ i wyedytowaæ plik setup.h ustawiaj¹c:
   - wxUSE_GLCANVAS na 1
   - wxUSE_GUI na 1
   - wxUSE_UNICODE 0
   - opcjonalnie wxXLOCALS na 0 (wxWidgets domyslnie u¿ywaj¹ ustawieñ regionalnych, np ','  zamiast '.')
4. uruchomiæ wiersz polecenia œrodowsika kompilatora Visual Studio (Visual Studio 2008 x64 Win64 Command Prompt).
5. wejœæ do katalogu sciezka\do\wxWidgets\build\msw i uruchomic polecenie nmake z parametrami:
   - wersja Debug: 
	nmake -f makefile.vc TARGET_CPU=AMD64 BUILD=debug UNICODE=0 SHARED=1 USE_OPENGL=1 USE_GUI=1
   - wersja Release: 
	nmake -f makefile.vc TARGET_CPU=AMD64 BUILD=release UNICODE=0 SHARED=1 USE_OPENGL=1 USE_GUI=1
6. po kompilacji sprawdziæ czy w katalogu œciezka\do\wxWidgets\lib\vc_amd64_dll znajduj¹ siê biblioteki:
   - wersja Debug: wxbase291d_vc_custom.dll wxmsw291d_core_vc_custom.dll wxmsw291d_gl_vc_custom.dll
   - wersja Debug: wxbase291_vc_custom.dll wxmsw291_core_vc_custom.dll wxmsw291_gl_vc_custom.dll
7. powy¿sze biblioteki skopiowac do odpowiednich katalogów uruchomieniowych.
8. w katalogach tych powinna sie znaleŸæ biblioteka glut32.dll (wersja x64).
9. po uruchomieniu programu mo¿emy bezpoœrednio wczytywac pliki, tzn. femViewer korzysta z wewnetrznych
   procedur (tylko siatki pryzmatyczne) b¹dŸ skorzystaæ z zewnêtrznych modu³ów jak modu³ siatek hybrydowych.
   W tym celu nale¿y przed wczytaniem pliku siatki wejœæ do menu Configuration i wybrac Use external approximation,
   co spowoduje ¿e domyslnym modu³em jest modu³ siatek hybrydowych.
10 dla siatek hybrydowych pliki field musz¹ zawierac w swojej nazwie cz³on fieldh_std.
11 dla siatek pryzmatycznych: dla DG field_dg, dla STD field_std.
10. pewne funkcje programu nie dzia³aj¹ i wkrótce zostanie naprawione.

Dla systemu Linux
proszê byc cierpliwym ;).   