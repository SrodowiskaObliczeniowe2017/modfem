FemViewr - aplikacja graficzna do wizualizacji rozwi�za� MES na siatkach pryzmatycznych i czworo�ciennych
Autor: Pawe� Macio�, Politechnika Krakowska
wxWidgtes-version: 2.9.1.
Kompilator: Microsoft VC++ 2008.
System testowy: Windows 7.

FemViewer dzia�a w oparciu o biblioteki wxWidgets, w tym celu nalezy mie� skompilowane biblioteki 
wxWidgets w postaci dynamicznej (.dll w Windows).

Przygotowanie (Win 7):
1. pobra� i zainstalowa� najnowsz� wersje bibliotek wxWidgets dla Win7.
2. w zmienych systemowych ustawic zmien� WXWIN na scie�k�\do\katalogu\wxWidgets (np. C:\wxWidgets-2.9.1)
3. wej�� do katalogu sciezka\do\wxWidgets\include\wx\msw\ i wyedytowa� plik setup.h ustawiaj�c:
   - wxUSE_GLCANVAS na 1
   - wxUSE_GUI na 1
   - wxUSE_UNICODE 0
   - opcjonalnie wxXLOCALS na 0 (wxWidgets domyslnie u�ywaj� ustawie� regionalnych, np ','  zamiast '.')
4. uruchomi� wiersz polecenia �rodowsika kompilatora Visual Studio (Visual Studio 2008 x64 Win64 Command Prompt).
5. wej�� do katalogu sciezka\do\wxWidgets\build\msw i uruchomic polecenie nmake z parametrami:
   - wersja Debug: 
	nmake -f makefile.vc TARGET_CPU=AMD64 BUILD=debug UNICODE=0 SHARED=1 USE_OPENGL=1 USE_GUI=1
   - wersja Release: 
	nmake -f makefile.vc TARGET_CPU=AMD64 BUILD=release UNICODE=0 SHARED=1 USE_OPENGL=1 USE_GUI=1
6. po kompilacji sprawdzi� czy w katalogu �ciezka\do\wxWidgets\lib\vc_amd64_dll znajduj� si� biblioteki:
   - wersja Debug: wxbase291d_vc_custom.dll wxmsw291d_core_vc_custom.dll wxmsw291d_gl_vc_custom.dll
   - wersja Debug: wxbase291_vc_custom.dll wxmsw291_core_vc_custom.dll wxmsw291_gl_vc_custom.dll
7. powy�sze biblioteki skopiowac do odpowiednich katalog�w uruchomieniowych.
8. w katalogach tych powinna sie znale�� biblioteka glut32.dll (wersja x64).
9. po uruchomieniu programu mo�emy bezpo�rednio wczytywac pliki, tzn. femViewer korzysta z wewnetrznych
   procedur (tylko siatki pryzmatyczne) b�d� skorzysta� z zewn�trznych modu��w jak modu� siatek hybrydowych.
   W tym celu nale�y przed wczytaniem pliku siatki wej�� do menu Configuration i wybrac Use external approximation,
   co spowoduje �e domyslnym modu�em jest modu� siatek hybrydowych.
10 dla siatek hybrydowych pliki field musz� zawierac w swojej nazwie cz�on fieldh_std.
11 dla siatek pryzmatycznych: dla DG field_dg, dla STD field_std.
10. pewne funkcje programu nie dzia�aj� i wkr�tce zostanie naprawione.

Dla systemu Linux
prosz� byc cierpliwym ;).   