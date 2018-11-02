#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string>
#include <windows.h>
#include <iostream>

using namespace std;

#define sign(x) ( x < 0 ? -1 : ( x == 0 ? 0 : 1) )
#define abs(x) ( x < 0 ? -(x) : (x) )
//#define min(x,y) ((x)<(y)?(x):(y))
//#define max(x,y) ((x)>(y)?(x):(y))
#define PI 3.14159265358979
#define INF 10000000
#define Q 520															// rozmiar tablic (kontrola ilości pamięci)


double psi[Q][Q], omega[Q][Q];											// niewiadome pola
double v_x[Q][Q], v_y[Q][Q], v_modul[Q][Q];								// pola predkości (wynikowe)
double r[Q][Q], e[Q][Q];												// wektory błędu i normy do multigrida
double temp[Q][Q], k_next[Q][Q], k1[Q][Q], k2[Q][Q], k3[Q][Q], k4[Q][Q];// do runge-kutty dla omega
double x[Q][Q], y[Q][Q], xp[Q][Q], yp[Q][Q];							// położenia dla animacji
unsigned char img[54 + 3 * Q*Q];										// tablica do funkcji tworzącej .bmp

// DANE INICJACJI (MOZNA PODDAWAC ZMIANOM):
int n = 20;																// ilość podziałów boku opływanego kradratu o boku 1
double Re = 250;														// liczba Reynoldsa (policzona dla prędkości strumienia na wlocie
double dt = 0.001;														// (?? może być za duża) długość jednego kroku czasowego
double eps = 1e-4;                                                      // dokladnosc liczenia dyskretyzacji w przestrzeni
int N = 200 + 1;														// ilość pól siatki w kierunku poziomym
int M = 100 + 1;														// ilość pól siatki w kierunku pionowym
double skala = 10.;                                                     // im wieksza wartosc 'skala', tym mocniej mniejsza roznica wartosci zmienia kolor
                                                                        // ^ efekt widac w benchamarku 'laplasjan.bmp'


// INNE DANE (NIE ZMIENIAĆ):
double a = (double)((double)N / n);										// długość prostokąta (podziel ten bok na a * n fragmentów)
double b = (double)((double)M / n);										// szerokość prostokąta (podziel ten bok na b * n fragmentów)
double x_k = (a - 1) / 4.;
double y_k = (b - 1) / 2.;												// współrzędne lewego dolnego rogu opływanego kwadratu
int i_k = (int)(x_k * (double)n + 1);
int j_k = (int)(y_k * (double)n + 1);									// indeksy KRATKI odpowiadającej lewemu dolnemu rogowi kwadratu
double h = (double)1 / n;												// krok odległościowy
double ni = 0.001;														// lepkość kinematyczna
double U = Re * ni;														// prędkość górnej ścianki
double T_max = 80. / U;													// czas charakterystyczny zagadnienia (na podstawie równości liczb Strouhala)




void hsv2rgb(double hue, double sat, double val, double &red, double &grn, double &blu) {   // zamiana koloru HSV na kolor RGB

	double i, f, p, q, t;

	if (val == 0) {
		red = 0;
		grn = 0;
		blu = 0;
	}
	else {
		hue /= 60;
		i = floor(hue);
		f = hue - i;
		p = val*(1 - sat);
		q = val*(1 - (sat*f));
		t = val*(1 - (sat*(1 - f)));
		if (i == 0) { red = val; grn = t; blu = p; }
		else if (i == 1) { red = q; grn = val; blu = p; }
		else if (i == 2) { red = p; grn = val; blu = t; }
		else if (i == 3) { red = p; grn = q; blu = val; }
		else if (i == 4) { red = t; grn = p; blu = val; }
		else if (i == 5) { red = val; grn = p; blu = q; }
	}
}


void wartosci_poczatkowe() {				// początkowe wartości (pierwsze przybliżenie) pól psi(x,y) i omega(x,y)

	int i, j;

	// warunki początkowe psi w obszarze i brzegowe psi na brzegu DOLNYM, i poczatkowe omega w obszarze i na brzegu:

	for (j = 1; j <= M ; j++)
		for (i = 1; i <= N; i++)
			psi[i][j] = omega[i][j] = 0;


	// warunki brzegowe psi na brzegu GÓRNYM:

	for (i = 1; i <= N; i++)
		psi[i][M] = U * b;


	// warunki brzegowe psi na wlocie i wylocie:

	for (j = 1; j <= M; j++)
		psi[1][j] = psi[N][j] = U * b * (double)(j-1) / (M-1);



	// wyzeruj wartosci WEWNĄTRZ kwadratu:
	for (j = j_k + 1; j < j_k + n - 1; j++)
		for (i = i_k + 1; i < i_k + n - 1; i++)
			;//psi[i][j] = 0;


	// zdefiniuj funkcję prądu na brzegach kwadratu (wynika z fizycznego sensu funkcji prądu!)
	for (i = 0; i < n; i++) {
		psi[i_k][j_k + i]			= U*b/2;
		psi[i_k + n - 1][j_k + i]	= U*b/2;
		psi[i_k + i][j_k]			= U*b/2;
		psi[i_k + i][j_k + n - 1]	= U*b/2;
	}
}

void wartosci_poczatkowe_z_pliku( string s) {	// wczytaj dane z poprzednio utworzonego pliku

	int i, j;
	FILE *p = fopen( s.c_str(), "r");

	// wczytaj początkowe psi:
	for (j = 1; j <= N;j++)
		for (i = 1; i <= N; i++)
			fscanf(p, "%lf", &psi[i][j]);

	// wczytaj początkowe omega:
	for (j = 1; j <= N;j++)
		for (i = 1; i <= N; i++)
			fscanf(p, "%lf", &omega[i][j]);

	fclose(p);
}

void wartosci_poczatkowe_z_pliku_i_rozszerzenie_siatki( string s) { // wczytaj siatke rzadszą i interpolują ją na gęstszą

	int i, j;
	FILE *p = fopen(s.c_str(), "r");

	// wczytaj wartości odpowiadające:
	for (j = 1; j <= N/2+1; j++)
		for (i = 1; i <= N/2+1; i++)
			fscanf(p, "%lf", &psi[2*i-1][2*j-1]);

	// uzupełnij wartości pośrednie w rzędach: (dla każdej nieparzystej kolumny, co drugi rząd trzeba uzupełnić)
	for (j = 1; j <= N; j+=2)
		for (i = 2; i <= N; i+=2)
			psi[i][j] = ( psi[i-1][j] + psi[i+1][j]) / 2;

	// uzupełnij wartości pośrednie w kolumnach: (dla każdego parzystego rzędu, KAŻDY element trzeba uzupełnić)
	for (j = 2; j <= N; j += 2)
		for (i = 1; i <= N; i++)
			psi[i][j] = (psi[i][j-1] + psi[i][j+1]) / 2;



	// wczytaj wartości odpowiadające:
	for (j = 1; j <= N / 2 + 1; j++)
		for (i = 1; i <= N / 2 + 1; i++)
			fscanf(p, "%lf", &omega[2 * i - 1][2 * j - 1]);

	// uzupełnij wartości pośrednie w rzędach: (dla każdej nieparzystej kolumny, co drugi rząd trzeba uzupełnić)
	for (j = 1; j <= N; j += 2)
		for (i = 2; i <= N; i += 2)
			omega[i][j] = (omega[i - 1][j] + omega[i + 1][j]) / 2;

	// uzupełnij wartości pośrednie w kolumnach: (dla każdego parzystego rzędu, KAŻDY element trzeba uzupełnić)
	for (j = 2; j <= N; j += 2)
		for (i = 1; i <= N; i++)
			omega[i][j] = (omega[i][j - 1] + omega[i][j + 1]) / 2;





	fclose(p);
}

void wypisz_do_pliku(string s) {

	int i, j;
	FILE *f = fopen(s.c_str(), "w");


	// wypisz psi do pliku:
	for (j = 1; j <= N; j++) {

		for (i = 1; i <= N; i++)
			fprintf(f, "%lf ", psi[i][j]);

		putc(10, f);
	}

	// wypisz omega do pliku:
	for (j = 1; j <= N; j++) {

		for (i = 1; i <= N; i++)
			fprintf(f, "%lf ", omega[i][j]);

		putc(10, f);
	}
	fclose(f);
}

bool czy_lezy_poza_kwadratem( int i, int j) {

	if (i_k <= i && i <= i_k + n - 1 && j_k <= j && j <= j_k + n - 1)
	//if (( i_k <= i && i <= i_k + n - 1 && ( j_k == j || j == j_k + n - 1) ) || ( (i_k == i || i == i_k + n - 1) && j_k <= j && j <= j_k + n - 1) )
		return 0;
	else
		return 1;
}

void rysuj_kolorowy_wykres(string s, double pole[Q][Q]) {		// kolorowy wykres danego pola

	int i, j;
	double mn = INF, mx = -INF, MN = INF;



	int w, h, x, y;
	w = N;
	h = M;

	FILE *f;

	int filesize = 54 + 3 * w*h;  //w is your image width, h is image height, both int




								  // znajdz ekstremalne wartosci pola
	for (i = 1; i <= N; i++)
		for (j = 1; j <= M; j++)
			mn = min(mn, pole[i][j]),
			mx = max(mx, pole[i][j]);

	printf("%min = %lf\nmax = %lf\n\n", mn, mx);           // wypisz ekstremalne wartosci rysowanego pola
	mn /= skala;                                           // modyfikujac te wartosc zmieniamy skale kolorow
	mx /= skala;                                           // modyfikujac te wartosc zmieniamy skale kolorow


	// czerwony kolor - maksymalna wartość
	// fioletowy kolor - minimalna wartość

	double r, g, b, hue;

	for (j = 1; j <= M; j++)
		for (i = 1; i <= N; i++) {

			hue = pole[i][j];
			hue = max(hue, mn);
			hue = min(hue, mx);

			hue = 360 - 360 * (hue - mn) / (mx - mn);
			hue *= 0.8;
			//hue = 360 - 360 * (hue - mn) / (mx - mn);

			while (hue > 360)
				hue -= 360;
			while (hue < 0)
				hue += 360;


			if(!czy_lezy_poza_kwadratem(i,j))
				hsv2rgb(hue, 1, 0, r, g, b);
			else
				hsv2rgb(hue, 1, 1, r, g, b);

			r *= 255;
			g *= 255;
			b *= 255;

			//filledrectangle(dx*(i - 1), dx*(j - 1), dx*i, dx*j);

			x = i - 1; y = (h - 1) - j + 1;
			img[(x + y*w) * 3 + 2] = (unsigned char)(r);
			img[(x + y*w) * 3 + 1] = (unsigned char)(g);
			img[(x + y*w) * 3 + 0] = (unsigned char)(b);
		}





	unsigned char bmpfileheader[14] = { 'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0 };
	unsigned char bmpinfoheader[40] = { 40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0 };
	unsigned char bmppad[3] = { 0,0,0 };

	bmpfileheader[2] = (unsigned char)(filesize);
	bmpfileheader[3] = (unsigned char)(filesize >> 8);
	bmpfileheader[4] = (unsigned char)(filesize >> 16);
	bmpfileheader[5] = (unsigned char)(filesize >> 24);

	bmpinfoheader[4] = (unsigned char)(w);
	bmpinfoheader[5] = (unsigned char)(w >> 8);
	bmpinfoheader[6] = (unsigned char)(w >> 16);
	bmpinfoheader[7] = (unsigned char)(w >> 24);
	bmpinfoheader[8] = (unsigned char)(h);
	bmpinfoheader[9] = (unsigned char)(h >> 8);
	bmpinfoheader[10] = (unsigned char)(h >> 16);
	bmpinfoheader[11] = (unsigned char)(h >> 24);

	f = fopen(s.c_str(), "wb");
	fwrite(bmpfileheader, 1, 14, f);
	fwrite(bmpinfoheader, 1, 40, f);
	for (i = 0; i < h; i++)
	{
		fwrite(img + (w*(h - i - 1) * 3), 3, w, f);
		fwrite(bmppad, 1, (4 - (w * 3) % 4) % 4, f);
	}
	fclose(f);
}

void stablicuj_pola_predkosci() {

	int i, j;

	// prędkości na brzegach górnym i dolnym obszaru (0, bo lepkość!)
	for (i = 1; i <= N; i++)
		v_x[i][1] = 0,
		v_y[i][1] = 0,
		v_x[i][M] = 0,
		v_y[i][M] = 0;


	// predkości na wlocie i wylocie:
	for (j = 2; j <= M-1; j++)
		v_x[1][j] = U,
		v_y[1][j] = 0,
		v_x[N][j] = U,
		v_y[N][j] = 0;




	for (j = 2; j <= M - 1; j++)
		for (i = 2; i <= N - 1; i++)
			if (czy_lezy_poza_kwadratem(i, j))
				v_x[i][j] = (psi[i][j + 1] - psi[i][j - 1]) / (2 * h),
				v_y[i][j] = -(psi[i + 1][j] - psi[i - 1][j]) / (2 * h);
			else
				v_x[i][j] = 0,
				v_y[i][j] = 0;


	for (j = 1; j <= M; j++)
		for (i = 1; i <= N; i++)
				v_modul[i][j] = sqrt(v_x[i][j] * v_x[i][j] + v_y[i][j] * v_y[i][j]);
}


double wygladz( double ilosc_iteracji, double eps) {	// po operacji na rzadkiej siatce multigrida, popraw wartości w kratkach interpolowanych


	int i, j, k;
	double u_p, alfa = -0.5, psi_old, sum_k;

	sum_k = INF;
	// wykonaj 3 iteracje G-S
	for (k = 0; sqrt(sum_k) > eps && k < ilosc_iteracji; k++) {

		sum_k = 0.;
		// dla wszystkich punktow wewnetrznych obszaru wykonaj: ... + RELAKSACJA
		for (j = 2; j <= M - 1; j++)
			for (i = 2; i <= N - 1; i++)
			if(czy_lezy_poza_kwadratem(i,j)){

				psi_old = psi[i][j];

				u_p = (psi[i - 1][j] + psi[i + 1][j] + psi[i][j - 1] + psi[i][j + 1] + h*h * omega[i][j]) / 4;
				psi[i][j] = alfa * psi[i][j] + (1 - alfa) * u_p;

				sum_k += (psi[i][j] - psi_old)*(psi[i][j] - psi_old);
			}
	}

	return sqrt(sum_k);
}

double multigrid( int m, int ilosc_iteracji, double eps) {

	// popraw dokładność rozwiązania: lap(psi) = - omega ; zwróć różnicę nowego psi w stosunku do psi z poprzedniego kroku ;

											// m - 2^stopień multigrida ( m musi być potęgą 2 )

											// Uwaga0: siatka rzadsza ma TEN SAM BRZEG co siatka wyjściowa

											// Uwaga1: pole psi nie wymaga korekcji na brzegu! zatem e_k[na brzegu] = 0. !

											// Uwaga2: ponieważ A*e_k = r_k, zaś e_k[na brzegu] = 0, to r_k[na brzegu] = 0. !

	int i, j, k, p;
	double sum_k, e_old;


	// znajdź wektor residuów r na siatce pełnej

	for (i = 2; i <= N - 1; i++)
		for (j = 2; j <= M - 1; j++)
		if ( czy_lezy_poza_kwadratem(i,j)) {

			temp[i][j] = r[i][j] = -omega[i][j] - (psi[i + 1][j] + psi[i - 1][j] + psi[i][j + 1] + psi[i][j - 1] - 4 * psi[i][j]) / (h*h);
			e[i][j] = 0.;		// inicjacja: wektor błędu wewnątrz obszaru = 0.
		}
		else {

			temp[i][j] = r[i][j] = e[i][j] = 0.;
		}




	// stwórz wnętrze rzadszej siatki metodą "wagową" (brzeg nie jest potrzebny)
	// jeśli poziom multigrida > 1, to zastosuj rekurencyjne przechodzenie do coraz niższych poziomów:
	for( p = 2 ; p <= m; p *= 2)
		for (i = 1+p; i <= N - p; i += p)
			for (j = 1+p; j <= M - p; j += p)
				temp[i][j] =
				(r[i - p/2][j - p/2] + r[i - p/2][j + p/2] + r[i + p/2][j - p/2] + r[i + p/2][j + p/2] + 2 * (r[i][j - p/2] + r[i][j + p/2] + r[i - p/2][j] + r[i + p/2][j]) + 4 * r[i][j]) / 16;




	// oblicz wektor błędu e[][] siatki rzadkiej z dokładnością do 1e-5
	k = 0;
	for (sum_k = INF; sqrt(sum_k) > eps || k < 3; k++) {

		sum_k = 0;

		for (i = 1 + m; i <= N - m; i += m)
			for (j = 1 + m; j <= M - m; j += m)
			if (czy_lezy_poza_kwadratem(i,j)){

				e_old = e[i][j];
				e[i][j] = (e[i - m][j] + e[i + m][j] + e[i][j - m] + e[i][j + m] - m*m * h*h * temp[i][j]) / 4;
				sum_k += (e[i][j] - e_old)*(e[i][j] - e_old);
			}
	}

	// ilość operacji do uzyskania zbieżności:
	printf("%d\n", k);



	// interpoluj e[][] (uzupełnij brakujące wartości, nie zmieniaj brzegu [dlatego start od i,j = 1] )
	// jeśli poziom multigrida > 1, to zastosuj rekurencyjne przechodzenie do coraz NIŻSZYCH poziomów:
	for (p = m; p >= 2; p /= 2)
		for (i = 1; i <= N - p; i += p)
			for (j = 1; j <= M - p; j += p) {

				e[i + p/2][j] = (e[i][j] + e[i + p][j]) / 2;
				e[i][j + p/2] = (e[i][j] + e[i][j + p]) / 2;
				e[i + p/2][j + p/2] = (e[i][j] + e[i + p][j] + e[i][j + p] + e[i + p][j + p]) / 4;
			}





	// popraw pole psi o obliczony błąd! Oblicz i zwróć błąd powstałego psi w stosunku do psi wyjœciowego
	// nie ma potrzeby poprawiania brzegu (dlatego start od i,j = 2)

	for (i = 2; i <= N - 1; i++)
		for (j = 2; j <= M - 1; j++)
			if(czy_lezy_poza_kwadratem(i,j))
			psi[i][j] += e[i][j];



	// A TERAZ policz błąd na samych punktach rzadkiej siatki:
	sum_k = 0.;
	for (i = 1+m; i <= N - m; i+=m)
		for (j = 1+m; j <= M - m; j+=m)
		if(czy_lezy_poza_kwadratem(i,j)){

			r[i][j] = -omega[i][j] - (psi[i + m][j] + psi[i - m][j] + psi[i][j + m] + psi[i][j - m] - 4 * psi[i][j]) / (h*h);
			sum_k += r[i][j] * r[i][j];
		}


	// "ilosc_iteracji" razy użyj G-S na pełnej siatce
	wygladz( ilosc_iteracji, eps);

	// zwróć błąd rzadkiej siatki:
	return sqrt(sum_k);
}

void rozwiaz_laplasjan_psi(double alfa, double eps) {	// oblicz pole psi z równania: laplasjan psi = - omega, alfa - wsp. relaksacji

											// dokładność do eps

	int i, j;
	double sum_k;



	// zadaj warunki brzegowe i "zabrzegowe" (ekseperyment) polom r i e, przed użyciem multigrida:

	for (i = 0; i <= N; i++)
		for (j = 0; j <= M; j++)
			r[i][j] = e[i][j] = 0.;


	// N = 128 + 1
	//multigrid(8, 8, eps);
	//multigrid(4, 4, eps);
	//multigrid(2, 4, eps);
	//multigrid(1, 1, eps);



	// N = 256 + 1
	//multigrid(16, 10, eps);
	//multigrid(8, 8, eps);
	//multigrid(4, 6, eps);
	//multigrid(2, 4, eps);
	//multigrid(1, 1, eps);




	// N = 512 + 1
	//multigrid(32, 25, eps);
	//multigrid(16, 10, eps);
	//multigrid(8, 5, eps);
	//multigrid(4, 5, eps);
	//multigrid(2, 3, eps);
	//multigrid(1, 1, eps);

	multigrid(16, 10, eps);
	multigrid(8, 8, eps);
	multigrid(4, 6, eps);
	multigrid(2, 4, eps);
	multigrid(1, 1, eps);





	// policz błąd2 (odchylenie standardowe normy):
	sum_k = 0.;
	for (i = 2; i <= N - 1; i++)
		for (j = 2; j <= M - 1; j++) {

			r[i][j] = -omega[i][j] - (psi[i + 1][j] + psi[i - 1][j] + psi[i][j + 1] + psi[i][j - 1] - 4 * psi[i][j]) / (h*h);
			sum_k += r[i][j] * r[i][j];
		}


	//printf("blad1 = %lf ", wygladz(1));
	//printf("blad2 = %lf\n", sqrt(sum_k) / N);

}

void rozwiaz_laplasjan_psi2(double alfa, double eps) {	// BRUTEFORCE. DZIAŁA

	int i, j, k;
	double u_p, sum_k, psi_old;

	sum_k = INF;
	k = 0;
	// uzyskaj jakąś dokładność za pomocą G-S (1e-3)
	while (sum_k > eps) {


		sum_k = 0.;
		// dla wszystkich punktow wewnetrznych obszaru wykonaj: ... + RELAKSACJA

		for (j = 2; j <= M - 1; j++) {
			for (i = 2; i <= N - 1; i++)
				if( czy_lezy_poza_kwadratem(i,j)) {

					u_p = (psi[i - 1][j] + psi[i + 1][j] + psi[i][j - 1] + psi[i][j + 1] + h*h * omega[i][j]) / 4;

					psi_old = psi[i][j];

					psi[i][j] = alfa * psi[i][j] + (1 - alfa) * u_p;

					sum_k += (psi[i][j] - psi_old)*(psi[i][j] - psi_old);
				}
		}


		sum_k = sqrt(sum_k);	// dokładność obliczenia pola psi z równania: laplasjan psi = - omega
		k++;
	}



	//printf("%d\n", k);

}

void runge_kutta(double k_next[Q][Q], double k[Q][Q]) {				// oblicz:	k_next = h * f(a * k)

		int i, j;
		double u_p;


		for (i = 2; i <= N - 1; i++)
			for (j = 2; j <= M - 1; j++)
				if( czy_lezy_poza_kwadratem( i, j)) {

					u_p =
						-(psi[i][j + 1] - psi[i][j - 1])*(k[i + 1][j] - k[i - 1][j]) / 4 +
						(psi[i + 1][j] - psi[i - 1][j])*(k[i][j + 1] - k[i][j - 1]) / 4 +
						ni * (k[i + 1][j] + k[i - 1][j] + k[i][j + 1] + k[i][j - 1] - 4 * k[i][j]);

					k_next[i][j] = dt * u_p / (h*h);
				}
}

void Gauss_Seidel(double alfa, double pole[Q][Q], int ilosc_slajdow, double eps) {

	// główna część algorytmu: parametr relaksacji ; które pole wypisywane ;  ile klatek ma się wykonać w tym czasie ; dokładność liczenia: lap psi = - omega


	int i, j, k, k_max;
	double T, t = GetTickCount();
	char str[80];


	// iteracje ukladu rownań:






	k_max = (int)(round(T_max / dt));

	for (k = 0, T = 0.; k <= k_max; k++, T += dt) {



		// poprawianie wirowości (omegi) na brzegach górnym i dolnym obszaru:
		for (i = 1; i <= N; i++)
			omega[i][1] = 2 / (h*h) * (psi[i][1] - psi[i][2]),
			omega[i][M] = 2 / (h*h) * (psi[i][M] - psi[i][M - 1]);


		// poprawianie wirowości (omegi) na wlocie i wylocie:
		for (j = 2; j <= M-1; j++)
			omega[1][j] = 2 / (h*h) * (psi[1][j] - psi[2][j]),		// wlot
			omega[N][j] = 2 / (h*h) * (-psi[N][j] + psi[N - 1][j]);	// wylot





		// popraw wirowość na wszystkich brzegach kwadratu!
		for (i = 0; i < n; i++) {
			omega[i_k][j_k + i]			 = 2 / (h*h) * (psi[i_k][j_k + i] - psi[i_k - 1][j_k + i]);			// lewy brzeg kwadratu
			omega[i_k + n-1][j_k + i]	 = 2 / (h*h) * (psi[i_k + n-1][j_k + i] - psi[i_k + n][j_k + i]);	// prawy brzeg kwadratu
			omega[i_k + i][j_k]			 = 2 / (h*h) * (psi[i_k + i][j_k] - psi[i_k + i][j_k - 1]);			// dolny brzeg kwadratu
			omega[i_k + i][j_k + n-1]	 = 2 / (h*h) * (psi[i_k + i][j_k + n-1] - psi[i_k + i][j_k + n]);	// górny brzeg kwadratu
		}



		// nowe pole omega (x,y):

		// zagadnienie czasowe:

		runge_kutta(k1, omega);			// k1 = h * f(omega)

		for (i = 2; i <= N - 1; i++)
			for (j = 2; j <= M - 1; j++)
				if( czy_lezy_poza_kwadratem( i,j))
					omega[i][j] += k1[i][j];





		rozwiaz_laplasjan_psi2(alfa, 1e-4);	// rozwiązuje: laplasjan psi = - omega , z dokładnością do 1e-6


		//printf("%d\n", k);




		if (k % (k_max / ilosc_slajdow) == 0) {				// ilość_klatek kroków czasowych w ciągu T = T_max

			printf("k = %d, U = %.2lf, T = %lf\n", k, U, T);    // k - numer kroku czasowego ; U - predkosc plynu na inlecie, T - aktualny czas

			sprintf(str, "animacja\\t = %.3lf.bmp", T);



			printf("czas kroku = %.3lfs\n", (GetTickCount() - t) / 1000);

			t = GetTickCount();


			stablicuj_pola_predkosci();

			rysuj_kolorowy_wykres( str, pole);
		}
	}
}

void benchmark( double eps) {	// rozwiąż: laplasjan(psi) = -omega = -(-1) = 1 na danej siatce N x M

	int i, j;


	// wypisz też dane przyszłej symulacji:

	printf("ilosc pol siatki w kierunku poziomym = %d\nilosc pol siatki w kierunku pionowym = %d\nRe = %.lf\n", N, M, Re);
	printf("dokladnosc wzgledem przestrzeni = %.10lf\ndokladnosc wzgledem czasu       = %.10lf\n\n", eps, dt);

	FILE *f = fopen("animacja\\dane.txt", "w");

	fprintf(f, "ilosc pol siatki w kierunku poziomym = %d\nilosc pol siatki w kierunku pionowym = %d\nRe = %.lf\n", N, M, Re);
	fprintf(f, "dokladnosc wzgledem przestrzeni = %.10lf\ndokladnosc wzgledem czasu       = %.10lf\n\n", eps, dt);

	fclose(f);




	// wartości początkowe psi zerowe

	for (j = 1; j <= N ; j++)
		for (i = 1; i <= M; i++)
			psi[i][j] = 0;

	// omega(x,y) = -1

	for (i = 1; i <= N; i++)
		for (j = 1; j <= M; j++)
			omega[i][j] = -1;


	// oblicz czas obliczenia równania laplasjanowego:

	double t = GetTickCount();



	rozwiaz_laplasjan_psi2(-0.5, eps);



	printf("czas benchmarku = %.3lfs\n", (GetTickCount() - t) / 1000);

	rysuj_kolorowy_wykres("animacja\\laplasjan.bmp", psi);
}

int main()
{
	// kwadrat o bokach: a (wzdłuż osi x) i b (wzdłuż osi y), z dolnym lewym wierzchołkiem w (0,0)

	int i, j;
	//double eps = 1e-4;

    CreateDirectory(TEXT("animacja"), NULL);



	// benchmark w postaci równania: laplasjan psi = 1
	benchmark(eps);
	//printf("a = %lf b = %lf\n", a, b);


	wartosci_poczatkowe();

	//printf("%d %d %lf %lf\n", i_k, j_k, x_k, y_k);





	Gauss_Seidel(-0.5, omega, 200, eps);        // 'omega' to pole, ktore bedzie wypisywane do pliku




	//wypisz_do_pliku( "xyz.txt");



	return 0;
}
