The program computes and draws an animation of vorticity of a flow around a square through a closed water channel at 50 < Re < 300. Von Karman vortex street is present.


=================================================================


Program oblicza przep�yw r�wnoleg�ego na wej�ciu strumienia p�ynu, przez prostok�tny kana� o sko�czonej d�ugo�ci w liczbie wymiar�w 2D. Profil pr�dko�ci na wylocie za�o�ony r�wnoleg�y. W odleg�o�ci po�owy szeroko�ci kana�u od jego wlotu znajduje si� przeszkoda w postaci kwadratu.

Po otworzeniu pliku main mamy mo�liwo�� modyfikacji ustawie� wej�ciowych, poprzez zmian� nast�puj�cych parametr�w globalnych:


=================================================================


int n = 20	/* Ilo�� podzia��w boku op�ywanego kradratu o boku 1.
Zwi�kszaj�c t� warto�� zwi�kszamy rozmiar kwadratu, przy zachowaniu tych samych wymiar�w kana�u */

double Re = 250;		/* Liczba Reynoldsa (policzona dla pr�dko�ci strumienia na wlocie */

double dt = 0.001;	/* D�ugo�� jednego kroku czasowego. Dla wi�kszych Re nale�y t� warto�� zmniejszy�, aby zachowa� zbie�no�� programu. */

double eps = 1e-4;    /* dokladnosc liczenia dyskretyzacji w przestrzeni. Zmniejszenie tej warto�ci mo�e zmniejszy� b��dy numeryczne. */
                                             
int N = 200 + 1;		// ilo�� p�l siatki w kierunku poziomym								
int M = 100 + 1;		// ilo�� p�l siatki w kierunku pionowym	

double skala = 10.;   /* im wieksza wartosc 'skala', tym mocniej mniejsza roznica wartosci p�l zmienia kolor w ko�cowych obrazkach. */
      
                                                                  
=================================================================


Program sam wylicza pr�dko�� strumienia wej�ciowego, lepko�� kinematyczn� i liczb� Stroughala zapewniaj�ce optymalne wyniki.	

Wymiary kana�u:

D�ugo�� boku kwadratu = 1 jednostka [j]
D�ugo�� kana�u 		= N / n [j]
Szeroko�� kana�u 	= M / n [j]



Wywo�anie komendy "benchmark(eps);" na pocz�tku programu powoduje rozpocz�cie oblicze� nast�puj�cego r�wnania r�niczkowego:

laplasjan(psi) = 1 na danej siatce N x M

Jako wynik zwraca czas oblicze� (do okna programu) oraz plik obrazkowy "laplasjan.png" (do folderu "animacja" utworzonego w folderze z plikem "main.cpp") przedstawiaj�cy pokolorowany obrazek pola "psi" (skala kolor�w zale�y od warto�ci parametru "skala").


=================================================================


W oknie programu podczas oblicze� wy�wietlane s� nast�puj�ce dane, zawsze w momencie zapisu klatki czasowej do pliku graficznego ".bmp":

k = 1600, U = 0.25, T = 1.600000
czas kroku = 4.406s
min = -22.230503
max = 18.383575


- k oznacza ilo�� obliczonych do tej pory krok�w czasowych
- U oznacza pr�dko�� wlotow� do kana�u w [j / s]
- T oznacza aktualny czas ( k * dt) zapisywanej klatki czasowej.
- min/max oznaczaj� minimaln�/maksymaln� warto�� danego pola w obszarze kana�u.

Klatki czasowe zapisywane s� w formacie "t = T.bmp" w folderze "animacja", automatycznie tworzonym w folderze z plikiem "main.cpp".

Standardowo zapisywane do plik�w graficznych s� pola wirowo�ci (omega). Mo�liwe jest zapisywanie r�wnie� p�l: 'psi' (funkcja pr�du) oraz sk�adowych pr�dko�ci ('v_x', 'v_y', 'v_modul'). W tym celu nale�y podmieni� nast�puj�c� komend� w bloku int main():

Gauss_Seidel(-0.5, omega, 200, eps);

na np.:

Gauss_Seidel(-0.5, v_x, 200, eps);
