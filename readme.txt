The program computes and draws an animation of vorticity of a flow around a square through a closed water channel at 50 < Re < 300. Von Karman vortex street is present.

The attached pdf file, although corcerns a little different problem, describes most of CFD algorithms implemented in this program.


=================================================================


Program oblicza przepływ równoległego na wejściu strumienia płynu, przez prostokątny kanał o skończonej długości w liczbie wymiarów 2D. Profil prędkości na wylocie założony równoległy. W odległości połowy szerokości kanału od jego wlotu znajduje się przeszkoda w postaci kwadratu.

Dołączony plik pdf, choć dotyczy nieco innego zagadnienia, zawiera szczegółowy opis zaimplementowanych w programie algorytmów.


=================================================================


Po otworzeniu pliku main mamy możliwość modyfikacji ustawień wejściowych, poprzez zmianę następujących parametrów globalnych:


int n = 20	/* Ilość podziałów boku opływanego kradratu o boku 1.
Zwiększając tę wartość zwiększamy rozmiar kwadratu, przy zachowaniu tych samych wymiarów kanału */

double Re = 250;		/* Liczba Reynoldsa (policzona dla prędkości strumienia na wlocie */

double dt = 0.001;	/* Długość jednego kroku czasowego. Dla większych Re należy tę wartość zmniejszyć, aby zachować zbieżność programu. */

double eps = 1e-4;    /* dokladnosc liczenia dyskretyzacji w przestrzeni. Zmniejszenie tej wartości może zmniejszyć błędy numeryczne. */
                                             
int N = 200 + 1;		// ilość pól siatki w kierunku poziomym								
int M = 100 + 1;		// ilość pól siatki w kierunku pionowym	

double skala = 10.;   /* im wieksza wartosc 'skala', tym mocniej mniejsza roznica wartosci pól zmienia kolor w końcowych obrazkach. */
      
                                                                  
=================================================================


Program sam wylicza prędkość strumienia wejściowego, lepkość kinematyczną i liczbę Stroughala zapewniające optymalne wyniki.	

Wymiary kanału:

Długość boku kwadratu = 1 jednostka [j]
Długość kanału 		= N / n [j]
Szerokość kanału 	= M / n [j]



Wywołanie komendy "benchmark(eps);" na początku programu powoduje rozpoczęcie obliczeń następującego równania różniczkowego:

laplasjan(psi) = 1 na danej siatce N x M

Jako wynik zwraca czas obliczeń (do okna programu) oraz plik obrazkowy "laplasjan.png" (do folderu "animacja" utworzonego w folderze z plikem "main.cpp") przedstawiający pokolorowany obrazek pola "psi" (skala kolorów zależy od wartości parametru "skala").


=================================================================


W oknie programu podczas obliczeń wyświetlane są następujące dane, zawsze w momencie zapisu klatki czasowej do pliku graficznego ".bmp":

k = 1600, U = 0.25, T = 1.600000
czas kroku = 4.406s
min = -22.230503
max = 18.383575


- k oznacza ilość obliczonych do tej pory kroków czasowych
- U oznacza prędkość wlotową do kanału w [j / s]
- T oznacza aktualny czas ( k * dt) zapisywanej klatki czasowej.
- min/max oznaczają minimalną/maksymalną wartość danego pola w obszarze kanału.

Klatki czasowe zapisywane są w formacie "t = T.bmp" w folderze "animacja", automatycznie tworzonym w folderze z plikiem "main.cpp".

Standardowo zapisywane do plików graficznych są pola wirowości (omega). Możliwe jest zapisywanie również pól: 'psi' (funkcja prądu) oraz składowych prędkości ('v_x', 'v_y', 'v_modul'). W tym celu należy podmienić następującą komendę w bloku int main():

Gauss_Seidel(-0.5, omega, 200, eps);

na np.:

Gauss_Seidel(-0.5, v_x, 200, eps);
