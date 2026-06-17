# numlib - Biblioteka Metod Numerycznych

Biblioteka **numlib** implementuje w języku C++ (standard C++17) podstawowe algorytmy numeryczne poznane podczas laboratorium z Metod Numerycznych. Wszystkie funkcje zebrane są w przestrzeni nazw `numlib`.

## Struktura projektu

```
numlib/
├── include/              # Pliki nagłówkowe (interfejs publiczny)
│   ├── numlib.h          # Główny plik - wystarczy dołączyć ten jeden
│   ├── linear_systems.h  # Rozwiązywanie układów równań liniowych (LU)
│   ├── interpolation.h   # Interpolacja (Newton, Lagrange, Horner)
│   ├── approximation.h   # Aproksymacja średniokwadratowa
│   ├── integration.h     # Całkowanie numeryczne
│   ├── ode.h             # Równania różniczkowe (Euler, Heun, RK4, ...)
│   └── nonlinear.h       # Równania nieliniowe (Newton, sieczne)
├── src/                  # Implementacje
│   ├── linear_systems.cpp
│   ├── interpolation.cpp
│   ├── approximation.cpp
│   ├── integration.cpp
│   ├── ode.cpp
│   └── nonlinear.cpp
├── tests/                # Testy jednostkowe
│   ├── test_linear.cpp        (5 testów)
│   ├── test_interpolation.cpp (6 testów)
│   ├── test_integration.cpp   (7 testów)
│   ├── test_ode.cpp           (7 testów)
│   ├── test_nonlinear.cpp     (7 testów)
│   └── test_approximation.cpp (5 testów)
├── examples/
│   ├── example_integration_nonlinear.cpp
│   └── example_ode_interpolation.cpp
├── CMakeLists.txt
└── README.md
```

## Wymagania

- Kompilator C++17: `g++ >= 7` lub `clang++ >= 5`
- CMake >= 3.10
- Biblioteka standardowa + `libm`

## Kompilacja i uruchomienie

### Sposób 1 - CMake

```bash
mkdir build && cd build
cmake ..
cmake --build .

# Uruchom testy przez CTest
ctest --output-on-failure

# Uruchom przykłady
./example1
./example2
```

### Sposób 2 - kompilacja ręczna

```bash
# Kompilacja biblioteki
g++ -std=c++17 -O2 -Iinclude -c src/*.cpp
ar rcs libnumlib.a *.o

# Własny program korzystający z biblioteki
g++ -std=c++17 -O2 -Iinclude moj_program.cpp libnumlib.a -lm -o moj_program
```
### Sposób 3 - CMake (Windows - Terminal PowerShell / CMD)

Wymagane jest posiadanie zainstalowanego kompilatora (np. MinGW-w64 lub MSVC) oraz programu CMake dodanego do zmiennych środowiskowych PATH.

```bash
mkdir build
cd build
cmake ..
cmake --build .

# Uruchom testy przez CTest
ctest --output-on-failure

# Uruchom przykłady
.\example1.exe
.\example2.exe
```
### Sposób 4 - Pełne Visual Studio (Windows - Natywne wsparcie CMake)

Metoda przeznaczona dla standardowego środowiska Microsoft Visual Studio (np. 2019 / 2022) z zainstalowanym pakietem "Programowanie aplikacji klasycznych w C++". Nie wymaga ręcznego wpisywania żadnych poleceń w terminalu.

#### 1. Otwarcie projektu
1. Uruchom Visual Studio.
2. W oknie startowym wybierz opcję **Otwórz folder** (ang. *Open a local folder*).
3. Wskaż główny katalog biblioteki (ten, w którym znajduje się plik `CMakeLists.txt`).
4. Opcjonalnie można otworzyć folder wybierając opcję **Otwórz w programie Visual Studio**

#### 2. Wybór celu i uruchomienie (Przykłady / Testy)
Kompilacja i uruchamianie programów odbywają się bezpośrednio z poziomu górnego paska narzędzi:

1. **Konfiguracja:** Po lewej stronie od przycisku uruchamiania upewnij się, że wybrany jest tryb architektury (np. `x64-Debug` lub `x64-Release`).
2. **Wybór programu (Cel startowy):** Kliknij na rozwijaną listę tuż obok zielonej strzałki (domyślnie: *Wybierz element startowy*) i wskaż cel, który chcesz skompilować:
   * Wybierz **`example1.exe`** lub **`example2.exe`**, aby uruchomić programy demonstracyjne.
   * Wybierz konkretny plik testowy (np. **`test_linear.exe`**), aby sprawdzić poprawność działania danego modułu.
3. **Uruchomienie:** Kliknij zieloną strzałkę **Lokalny debuger systemu Windows** (lub naciśnij klawisz `F5`). Program automatycznie skompiluje się, zlinkuje i otworzy okno konsoli z wynikami.

#### 3. Zbiorcze uruchamianie testów (Eksplorator testów)
Dzięki wbudowanej integracji systemu CTest z Visual Studio, wszystkimi testami można zarządzać graficznie:
1. Z górnego menu wybierz **Test -> Eksplorator testów** (ang. *Test Explorer*).
2. W panelu, który pojawi się po lewej stronie, wyświetli się pełna lista zaimplementowanych modułów testowych.
3. Kliknij ikonę **Uruchom wszystkie** (dwie zielone strzałki). Środowisko automatycznie zweryfikuje wszystkie asercje i podświetli poprawne testy na zielono.

## Przegląd funkcji

### 1. Rozwiązywanie układów równań liniowych (`linear_systems.h`)

| Funkcja | Opis |
|---|---|
| `lu_decompose(A, n, L, U)` | Rozkład A = L·U |
| `lu_solve(A, b, n)` | Rozwiązuje Ax = b metodą LU |

```cpp
#include "numlib.h"
// 2x + y = 5,  4x + 3y = 11  =>  x=2, y=1
std::vector<long double> A = {2,1, 4,3};
std::vector<long double> b = {5, 11};
auto x = numlib::lu_solve(A, b, 2);
// x[0] == 2.0,  x[1] == 1.0
```

### 2. Interpolacja (`interpolation.h`)

| Funkcja | Opis |
|---|---|
| `newton_coefficients(xi, fxi)` | Ilorazy różnicowe Newtona |
| `newton_evaluate(x, xi, coeffs)` | Wartość wielomianu Newtona |
| `lagrange_evaluate(x, xi, fxi)` | Wartość wielomianu Lagrange'a |
| `horner(a, x)` | Wartość wielomianu schematem Hornera |
| `polynomial_natural(a, x)` | Wartość wielomianu metodą naturalną |

```cpp
std::vector<long double> xi  = {0.0L, 1.0L, 2.0L};
std::vector<long double> fxi = {1.0L, 2.0L, 5.0L};  // f(x) = x^2 + 1
auto c = numlib::newton_coefficients(xi, fxi);
long double y = numlib::newton_evaluate(1.5L, xi, c); // 3.25
```

### 3. Aproksymacja (`approximation.h`)

| Funkcja | Opis |
|---|---|
| `polynomial_approximation(f, a, b, deg)` | Współczynniki wielomianu aproksymującego |
| `approx_evaluate(coeffs, x)` | Wartość wielomianu aproksymującego |

```cpp
auto c = numlib::polynomial_approximation(sinl, 0.0L, 3.14159L, 5);
long double y = numlib::approx_evaluate(c, 1.0L); // ≈ sin(1) = 0.8415
```

### 4. Całkowanie numeryczne (`integration.h`)

| Funkcja | Opis |
|---|---|
| `trapezoidal(f, a, b, n)` | Złożona metoda trapezów, błąd O(h²) |
| `simpson(f, a, b, n)` | Złożona metoda Simpsona 1/3, błąd O(h⁴) |
| `gauss_quadrature(f, a, b, nodes)` | Kwadratura Gaussa-Legendre'a (2, 3 lub 4 węzły) |
| `gauss_composite(f, a, b, nodes, m)` | Złożona kwadratura Gaussa (m podprzedziałów) |

```cpp
long double I = numlib::gauss_composite(f, -1.0L, 15.0L, 4, 50);
```

### 5. Równania różniczkowe (`ode.h`)

Wszystkie metody rozwiązują zagadnienie Cauchy'ego: **y'(t) = f(t, y),  y(t₀) = y₀**  
Zwracają `vector<pair<long double, long double>>` - pary (tᵢ, yᵢ).

| Funkcja | Rząd błędu | Opis |
|---|---|---|
| `euler(f, y0, t0, t_end, N)` | O(h) | Metoda Eulera |
| `heun(f, y0, t0, t_end, N)` | O(h²) | Metoda Heuna (ulepszony Euler) |
| `midpoint(f, y0, t0, t_end, N)` | O(h²) | Metoda punktu środkowego |
| `rk4(f, y0, t0, t_end, N)` | O(h⁴) | Runge-Kutta 4. rzędu |

```cpp
auto f   = [](long double t, long double y){ return -y; }; // y' = -y
auto sol = numlib::rk4(f, 1.0L, 0.0L, 5.0L, 1000);
// sol.back() = {5.0, e^{-5}}
```

### 6. Równania nieliniowe (`nonlinear.h`)

| Funkcja | Opis |
|---|---|
| `numerical_derivative(f, x)` | Pochodna numeryczna (różnice centralne) |
| `newton_analytic(f, df, x0)` | Metoda Newtona z analityczną pochodną |
| `newton_numeric(f, x0)` | Metoda Newtona z pochodną numeryczną |
| `secant(f, x0, x1)` | Metoda siecznych |

```cpp
auto f  = [](long double x){ return x*x - 2.0L; };
auto df = [](long double x){ return 2.0L*x; };
long double root = numlib::newton_analytic(f, df, 1.5L); // sqrt(2)
```

## Testowanie

Biblioteka zawiera łącznie **37 testów jednostkowych** (co najmniej 5 na każdą kategorię):

```
tests/test_linear.cpp        - 5 testów
tests/test_interpolation.cpp - 6 testów
tests/test_integration.cpp   - 7 testów
tests/test_ode.cpp           - 7 testów
tests/test_nonlinear.cpp     - 7 testów
tests/test_approximation.cpp - 5 testów
```

Testy obejmują zarówno przypadki poprawne (sprawdzenie wartości numerycznych), jak i przypadki błędne (macierz osobliwa -> pusty wynik, zerowa pochodna -> NaN).

## Autorzy

Projekt zrealizowany w ramach laboratorium z **Metod Numerycznych** - ćwiczenie 13.
