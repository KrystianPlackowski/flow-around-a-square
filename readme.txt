The program computes and draws an animation of vorticity of a flow around a square through a closed water channel at 50 < Re < 300. Von Karman vortex street is present.


=================================================================


Program oblicza przep³yw równoleg³ego na wejœciu strumienia p³ynu, przez prostok¹tny kana³ o skoñczonej d³ugoœci w liczbie wymiarów 2D. Profil prêdkoœci na wylocie za³o¿ony równoleg³y. W odleg³oœci po³owy szerokoœci kana³u od jego wlotu znajduje siê przeszkoda w postaci kwadratu.

Po otworzeniu pliku main mamy mo¿liwoœæ modyfikacji ustawieñ wejœciowych, poprzez zmianê nastêpuj¹cych parametrów globalnych:


=================================================================


int n = 20	/* Iloœæ podzia³ów boku op³ywanego kradratu o boku 1.
Zwiêkszaj¹c tê wartoœæ zwiêkszamy rozmiar kwadratu, przy zachowaniu tych samych wymiarów kana³u */

double Re = 250;		/* Liczba Reynoldsa (policzona dla prêdkoœci strumienia na wlocie */

double dt = 0.001;	/* D³ugoœæ jednego kroku czasowego. Dla wiêkszych Re nale¿y tê wartoœæ zmniejszyæ, aby zachowaæ zbie¿noœæ programu. */

double eps = 1e-4;    /* dokladnosc liczenia dyskretyzacji w przestrzeni. Zmniejszenie tej wartoœci mo¿e zmniejszyæ b³êdy numeryczne. */
                                             
int N = 200 + 1;		// iloœæ pól siatki w kierunku poziomym								
int M = 100 + 1;		// iloœæ pól siatki w kierunku pionowym	

double skala = 10.;   /* im wieksza wartosc 'skala', tym mocniej mniejsza roznica wartosci pól zmienia kolor w koñcowych obrazkach. */
      
                                                                  
=================================================================


Program sam wylicza prêdkoœæ strumienia wejœciowego, lepkoœæ kinematyczn¹ i liczbê Stroughala zapewniaj¹ce optymalne wyniki.	

Wymiary kana³u:

D³ugoœæ boku kwadratu = 1 jednostka [j]
D³ugoœæ kana³u 		= N / n [j]
Szerokoœæ kana³u 	= M / n [j]



Wywo³anie komendy "benchmark(eps);" na pocz¹tku programu powoduje rozpoczêcie obliczeñ nastêpuj¹cego równania ró¿niczkowego:

laplasjan(psi) = 1 na danej siatce N x M

Jako wynik zwraca czas obliczeñ (do okna programu) oraz plik obrazkowy "laplasjan.png" (do folderu "animacja" utworzonego w folderze z plikem "main.cpp") przedstawiaj¹cy pokolorowany obrazek pola "psi" (skala kolorów zale¿y od wartoœci parametru "skala").


=================================================================


W oknie programu podczas obliczeñ wyœwietlane s¹ nastêpuj¹ce dane, zawsze w momencie zapisu klatki czasowej do pliku graficznego ".bmp":

k = 1600, U = 0.25, T = 1.600000
czas kroku = 4.406s
min = -22.230503
max = 18.383575


- k oznacza iloœæ obliczonych do tej pory kroków czasowych
- U oznacza prêdkoœæ wlotow¹ do kana³u w [j / s]
- T oznacza aktualny czas ( k * dt) zapisywanej klatki czasowej.
- min/max oznaczaj¹ minimaln¹/maksymaln¹ wartoœæ danego pola w obszarze kana³u.

Klatki czasowe zapisywane s¹ w formacie "t = T.bmp" w folderze "animacja", automatycznie tworzonym w folderze z plikiem "main.cpp".

Standardowo zapisywane do plików graficznych s¹ pola wirowoœci (omega). Mo¿liwe jest zapisywanie równie¿ pól: 'psi' (funkcja pr¹du) oraz sk³adowych prêdkoœci ('v_x', 'v_y', 'v_modul'). W tym celu nale¿y podmieniæ nastêpuj¹c¹ komendê w bloku int main():

Gauss_Seidel(-0.5, omega, 200, eps);

na np.:

Gauss_Seidel(-0.5, v_x, 200, eps);
