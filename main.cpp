#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cmath>

double EPSILON = 10e-8;

class Matrica {
public:
    int rows, cols;
    double** elements;

public:
    /**
     * Konstruktor bez argumenata koji inicijalizira praznu matricu.
     */
    Matrica() : rows(0), cols(0), elements(nullptr) {}

    /**
     * Konstruktor koji učitava matricu iz datoteke.
     * @param filename Putanja do datoteke koja sadrži matricu.
     * Datoteka mora sadržavati matricu gdje su vrijednosti odvojene razmacima.
     */
    explicit Matrica(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Can't open the file.");
        }
        std::string line;
        std::vector<std::vector<double>> tmpMatrica;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::vector<double> row;
            double val;
            while (iss >> val) {
                row.push_back(val);
            }
            tmpMatrica.push_back(row);
        }
        rows = (int) tmpMatrica.size();
        cols = (int) tmpMatrica[0].size();
        elements = new double*[rows];
        for (int i = 0; i < rows; i++) {
            elements[i] = new double[cols];
            for (int j = 0; j < cols; ++j) {
                elements[i][j] = tmpMatrica[i][j];
            }
        }
        file.close();
    }

    /**
     * Konstruktor koji inicijalizira matricu zadanim brojem redaka i stupaca.
     * @param rows Broj redaka.
     * @param cols Broj stupaca.
     */
    Matrica(int rows, int cols) : rows(rows), cols(cols) {
        elements = new double*[rows];
        for (int i = 0; i < rows; ++i) {
            elements[i] = new double[cols];
        }
    }

    /**
     * Operator dodjele kopira elemente druge matrice u ovu matricu.
     * @param other Matrica koja se kopira.
     * @return Referenca na ovu matricu.
     */
    Matrica& operator=(const Matrica& other) {
        if (this == &other) {
            return *this;
        }
        for (int i = 0; i < rows; ++i) {
            delete[] elements[i];
        }
        delete[] elements;
        rows = other.rows;
        cols = other.cols;
        elements = new double*[rows];
        for (int i = 0; i < rows; ++i) {
            elements[i] = new double[cols];
            for (int j = 0; j < cols; ++j) {
                elements[i][j] = other.elements[i][j];
            }
        }
        return *this;
    }

    /**
     * Operator indeksiranja omogućuje pristup elementima matrice po indeksu.
     * @param index Indeks retka.
     * @return Niz koji predstavlja redak matrice.
     */
    double* operator[](int index) const {
        return elements[index];
    }

    /**
     * Operator usporedbe provjerava jesu li dvije matrice jednake unutar zadane tolerancije EPSILON.
     * @param other Matrica s kojom se uspoređuje.
     * @return True ako su matrice jednake, inače false.
     */
    bool operator==(const Matrica& other) const {
        if (rows != other.rows || cols != other.cols) {
            return false;
        }
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (std::abs(elements[i][j] - other.elements[i][j]) >= EPSILON) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Metoda za ispisivanje matrice na standardni izlaz.
     */
    void ispisi() const {
        for (int i = 0; i < rows; ++i) {
            std::cout << "[";
            for (int j = 0; j < cols; ++j) {
                std::cout << elements[i][j] << "\t";
            }
            std::cout << "]";
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    /**
     * Metoda koja transponira matricu.
     * @return Transponirana matrica.
     */
    [[nodiscard]] Matrica transponiraj() const {
        Matrica result(cols, rows);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.elements[j][i] = elements[i][j];
            }
        }
        return result;
    }

    /**
     * Metoda zamjenjuje redove u permutacijskoj matrici(vektoru).
     * Zamjenjuje vrijednosti između dva reda permutacijske matrice `P` na pozicijama `i` i `pivot`.
     * Koristi se tijekom postupka LUP dekompozicije za zamjenu redova matrice.
     *
     * @param P Permutacijska matrica u kojoj će se izvršiti zamjena.
     * @param i Indeks prvog retka koji će se zamijeniti.
     * @param pivot Indeks drugog retka koji će se zamijeniti.
     */
    static void zamijeni(Matrica &P, int i, int pivot) {
        if (i < 0 || i >= P.rows || pivot < 0 || pivot >= P.rows) {
            throw std::out_of_range("Indeksi izvan granica.");
        }
        std::swap(P[i][0], P[pivot][0]);
    }

    /**
     * Metoda za kopiranje matrice.
     * @return Kopija ove matrice.
     */
    [[nodiscard]] Matrica kopiraj() const {
        Matrica copy(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                copy.elements[i][j] = elements[i][j];
            }
        }
        return copy;
    }

    /**
     * Metoda za ispisivanje matrice u datoteku.
     * @param filename Putanja do datoteke u koju se ispisuje matrica.
     */
    void ispisiUDatoteku(const std::string& filename) const {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Ne mogu otvoriti datoteku za pisanje.");
        }
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                file << elements[i][j];
                if (j < cols - 1) {
                    file << " ";
                }
            }
            file << std::endl;
        }
        file.close();
    }

    /**
     * Operator za transponiranje matrice.
     * Omogućuje korištenje operatora '~' za dobivanje transponirane matrice.
     * @return Nova matrica koja je transponirana verzija ove matrice.
     *
     * Primjer:
     * Matrica A;
     * Matrica At = ~A;  // dobiva se transponirana matrica od A
     */
    Matrica operator~() const {
        return transponiraj();
    }

    /**
     * Operator za množenje matrice s realnim brojem.
     * @param scalar Skalar s kojim se množi matrica.
     * @return Nova matrica koja je rezultat množenja.
     */
    Matrica operator*(double scalar) const {
        if (std::fabs(scalar) < EPSILON) {
            throw std::runtime_error("Greška: Skalar je preblizu nuli.");
        }
        Matrica result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.elements[i][j] = elements[i][j] * scalar;
            }
        }
        return result;
    }

    /**
     * Operator množi dvije matrice.
     * @param other Matrica s kojom se množi (desni operand).
     * @return Rezultantna matrica nakon množenja.
     */
    Matrica operator*(const Matrica& other) const {
        if (cols != other.rows) {
            throw std::runtime_error("Greška pri množenju: Matrice nisu ulančane.");
        }
        Matrica result(rows, other.cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < other.cols; ++j) {
                result[i][j] = 0.;
                for (int k = 0; k < other.rows; ++k) {
                    result[i][j] += elements[i][k] * other.elements[k][j];
                }
            }
        }
        return result;
    }

    /**
     * Operator zbraja dvije matrice.
     * @param other Matrica koja se zbraja (desni operand).
     * @return Rezultantna matrica nakon zbrajanja.
     */
    Matrica operator+(const Matrica& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Greška pri zbrajanju: Matrice nisu jednakih dimenzija.");
        }
        Matrica result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result[i][j] = elements[i][j] + other.elements[i][j];
            }
        }
        return result;
    }

    /**
     * Operator oduzima dvije matrice.
     * @param other Matrica koja se oduzima (desni operand).
     * @return Rezultantna matrica nakon oduzimanja.
     */
    Matrica operator-(const Matrica& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Greška pri oduzimanju: Matrice nisu jednakih dimenzija.");
        }
        Matrica result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result[i][j] = elements[i][j] - other.elements[i][j];
            }
        }
        return result;
    }

    /**
     * Operator dodaje matricu drugoj matrici (zbrajanje u mjestu).
     * @param other Matrica koja se zbraja (desni operand).
     * @return Referenca na trenutnu matricu nakon zbrajanja.
     */
    Matrica& operator+=(const Matrica& other) {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Greška pri zbrajanju: Matrice nisu jednakih dimenzija.");
        }
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                elements[i][j] += other.elements[i][j];
            }
        }
        return *this;
    }

    /**
     * Operator oduzima matricu od druge matrice (oduzimanje u mjestu).
     * @param other Matrica koja se oduzima (desni operand).
     * @return Referenca na trenutnu matricu nakon oduzimanja.
     */
    Matrica& operator-=(const Matrica& other) {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Greška pri zbrajanju: Matrice nisu jednakih dimenzija.");
        }
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                elements[i][j] -= other.elements[i][j];
            }
        }
        return *this;
    }

    /**
     * Metoda mijenja dimenzije matrice na nove zadane dimenzije.
     * @param new_rows Broj redaka za novu matricu.
     * @param new_cols Broj stupaca za novu matricu.
     * @return Matrica s novim dimenzijama.
     */
    [[nodiscard]] Matrica promijeni_dimenzije(int new_rows, int new_cols) const{
        Matrica result(new_cols, new_cols);
        for (int i = 0; i < std::min(rows, new_rows); ++i) {
            for (int j = 0; j < std::min(cols, new_rows); ++j) {
                result.elements[i][j] = elements[i][j];
            }
        }
        return result;
    }

    /**
     * Metoda primjenjuje permutaciju na matricu prema zadanoj permutacijskoj matrici.
     * @param P Permutacijska matrica.
     * @return Nova matrica s primijenjenom permutacijom.
     */
    Matrica primijeni_permutaciju(Matrica &P) const {
        Matrica result(rows, cols);
        // Iteriranje po retcima i stupcima, permutira matricu A
        for (int i = 0; i < rows; ++i) {
            int new_index = (int)P[i][0];
            for (int j = 0; j < cols; ++j) {
                result.elements[i][j] = elements[new_index][j];
            }
        }
        return result;
    }

    /**
     * Metoda dohvaća određeni stupac matrice.
     * @param col_idx Indeks stupca koji treba dohvatiti.
     * @return Matrica koja predstavlja stupac kao vektor.
     */
    [[nodiscard]] Matrica dohvati_stupac(int col_idx) const {
        if (col_idx < 0 || col_idx >= cols) {
            throw std::out_of_range("Indeks stupca izvan granica.");
        }
        Matrica result(rows, 1);
        for (int i = 0; i < rows; ++i) {
            result[i][0] = elements[i][col_idx];
        }
        return result;
    }

    /**
     * Metoda postavlja vrijednosti stupca na temelju zadanog vektora.
     * @param col_idx Indeks stupca koji treba postaviti.
     * @param column Matrica koja sadrži vrijednosti za postavljanje stupca.
     */
    void postavi_stupac(int col_idx, const Matrica& column) const {
        // Provjera indeksa
        if (col_idx < 0 || col_idx >= cols) {
            throw std::out_of_range("Indeks stupca izvan granica.");
        }
        // Provjera dimenzija
        if (column.cols != 1 || column.rows != rows) {
            throw std::invalid_argument("Neispravne dimenzije vektora za postavljanje stupca.");
        }
        // Postavljanje
        for (int i = 0; i < rows; ++i) {
            elements[i][col_idx] = column.elements[i][0];
        }
    }

    /**
     * Metoda izvršava supstituciju unaprijed za rješavanje sustava linearnih jednadžbi s donjom trokutastom matricom.
     * @param b Matrica koja predstavlja vektor s desne strane jednadžbi.
     */
    void supstitucija_unaprijed(Matrica &b) const {
        // Provjera dimenzija - matrica mora biti kvadratna, a vektor s desne strane odgovarajućih dimenzija
        if (rows != cols || b.cols != 1 || b.rows != rows) {
            throw std::invalid_argument("Dimenzije matrice i vektora nisu kompatibilne za supstituciju unaprijed.");
        }
        int n = rows;

        // Ažuriranje vrijednosti u vektoru 'b' oduzimanjem već poznatih članova
        for (int i = 0; i <= n - 1; ++i) {
            // Ažuriranje vrijednosti u vektoru 'b' oduzimanjem već poznatih članova
            for (int j = i + 1; j < n; ++j) {
                b[j][0] -= elements[j][i] * b[i][0];
            }
        }
    }

    /**
     * Metoda izvršava supstituciju unatrag za rješavanje sustava linearnih jednadžbi s gornjetrokutastom matricom.
     * @param b Matrica koja predstavlja vektor s desne strane jednadžbi.
     *
     * @throws std::runtime_error Ako je pivot element preblizu nuli, što sprječava rješavanje.
     */
    void supstitucija_unatrag(Matrica &b) const {
        // Provjera dimenzija - matrica mora biti kvadratna, a vektor s desne strane odgovarajućih dimenzija
        if (rows != cols || b.cols != 1 || b.rows != rows) {
            throw std::invalid_argument("Dimenzije matrice i vektora nisu kompatibilne za supstituciju unaprijed.");
        }
        int n = rows;

        // Prolazak unatrag kroz jednadžbe
        for (int i = n - 1; i >= 0; --i) {
            // Provjera pivot elementa - ako je preblizu nuli, rješavanje nije moguće
            if (std::fabs(elements[i][i]) < EPSILON) {
                throw std::runtime_error("Pivot element je preblizu nuli tijekom supstitucije unatrag.");
            }

            // Podjela radi dobivanja konačne vrijednosti za tu jednadžbu
            b[i][0] /= elements[i][i];

            // Ažuriranje vrijednosti u vektoru 'b' za sve prethodne redove
            for (int j = 0; j < i; ++j) {
                b[j][0] -= elements[j][i] * b[i][0];
            }
        }
    }

    /**
     * Izvršava LU dekompoziciju matrice, dijeleći je na donjetrokutastu (L) i gornjetrokutastu (U) matricu.
     * @throws std::runtime_error Ako je pivot element preblizu nuli, što sprječava dekompoziciju.
     */
    void LU_dekompozicija() const {
        // Provjera je li matrica kvadratna
        if (rows != cols) {
            throw std::invalid_argument("LU dekompozicija zahtijeva kvadratnu matricu.");
        }
        int n = rows;

        // Iteracija kroz stupce za dekompoziciju
        for (int i = 0; i < n - 1; ++i) {
            // Provjera pivot elementa
            if (std::fabs(elements[i][i]) < EPSILON) {
                throw std::runtime_error("Pivot je preblizu nuli, LU dekompozicija nije moguća.");
            }

            // Izračun elemenata za donju trokutastu matricu (L)
            for (int j = i + 1; j < n; ++j) {
                elements[j][i] /= elements[i][i];

                // Ažuriranje vrijednosti u gornjetrokutastoj matrici (U)
                for (int k = i + 1; k < n; ++k) {
                    elements[j][k] -= elements[j][i] * elements[i][k];
                }
            }
        }
    }

    /**
     * Izvršava LUP dekompoziciju matrice, dodjeljujući permutacije i vraćajući broj zamjena redaka.
     * @param P Permutacijska matrica koja se popunjava.
     * @return Broj zamjena redaka tijekom dekompozicije.
     */
    int LUP_dekompozicija(Matrica &P) const {
        int n = rows;
        int swap = 0;

        // Provjera je li matrica kvadratna
        if (rows != cols) {
            throw std::invalid_argument("Matrica mora biti kvadratna za LUP dekompoziciju.");
        }

        // Inicijalizacija permutacijske matrice
        for (int i = 0; i < n; ++i) {
            P[i][0] = i;
        }

        int pivot;

        // Glavna petlja za dekompoziciju
        for (int i = 0; i < n - 1; ++i) {
            pivot = i;

            // Pronalaženje najvećeg elementa za pivotiranje
            for (int j = i + 1; j < n; ++j) {
                if (std::fabs(elements[(int)P[j][0]][i]) > std::fabs(elements[(int)P[pivot][0]][i])) {
                    pivot = j;
                }
            }

            // Zamjena redaka prema pivotiranju
            zamijeni(P, i, pivot);
            if (pivot != i) {
                swap++;
            }

            // Izračun elemenata za donju trokutastu matricu (L)
            for (int j = i + 1; j < n; ++j) {
                if(std::fabs(elements[(int)P[i][0]][i]) < EPSILON) {
                    continue;
                }

                elements[(int)P[j][0]][i] /= elements[(int)P[i][0]][i];

                // Ažuriranje vrijednosti u gornjetrokutastoj matrici (U)
                for (int k = i + 1; k < n; ++k) {
                    elements[(int)P[j][0]][k] -= elements[(int)P[j][0]][i] * elements[(int)P[i][0]][k];
                }
            }
        }
        return swap;
    }

    /**
     * Metoda za rješavanje sustava linearnih jednadžbi koristeći LU dekompoziciju.
     * @param b Vektor slobodnih članova.
     * @return Rješenje sustava.
     */
    Matrica rijesi_sustav_LU(Matrica &b) const {
        // Provjera kompatibilnosti dimenzija matrice i vektora
        if (rows != cols || b.cols != 1 || b.rows != rows) {
            throw std::invalid_argument("Dimenzije matrice i vektora nisu kompatibilne za rješavanje sustava.");
        }

        // Izvršavanje LU dekompozicije
        this->LU_dekompozicija();

        // Izvođenje supstitucija unaprijed i unatrag za rješavanje sustava
        this->supstitucija_unaprijed(b);
        this->supstitucija_unatrag(b);

        return b;
    }

    /**
     * Metoda rješava sustav linearnih jednadžbi koristeći LUP dekompoziciju.
     * @param b Matrica koja predstavlja vektor s desne strane jednadžbi.
     * @return Rješenje sustava linearnih jednadžbi kao vektor matrice.
     */
    Matrica rijesi_sustav_LUP(Matrica &b) {
        // Provjera kompatibilnosti dimenzija matrice i vektora
        if (rows != cols || b.cols != 1 || b.rows != rows) {
            throw std::invalid_argument("Dimenzije matrice i vektora nisu kompatibilne za rješavanje sustava.");
        }

        int n = rows;
        Matrica P(n, 1);

        // Izvršavanje LUP dekompozicije i permutacija
        LUP_dekompozicija(P);
        *this = primijeni_permutaciju(P);
        Matrica Pb = b.primijeni_permutaciju(P);

        // Izvođenje supstitucija unaprijed i unatrag za rješavanje sustava
        supstitucija_unaprijed(Pb);
        supstitucija_unatrag(Pb);

        return Pb;
    }

    /**
     * Metoda izračunava inverz matrice koristeći LUP dekompoziciju.
     * @return Inverzna matrica.
     * @throws std::invalid_argument Ako matrica nije kvadratna, inverzija nije moguća.
     */
    [[nodiscard]] Matrica inverz() const {
        // Provjera je li matrica kvadratna
        if (rows != cols) {
            throw std::invalid_argument("Matrica mora biti kvadratna za računanje inverzije LUP dekompozicijom.");
        }

        Matrica A_inv(rows, cols);
        int n = rows;
        Matrica A_cpy = this->kopiraj();
        Matrica P(n, 1);

        // Izvršavanje LUP dekompozicije i permutacija
        A_cpy.LUP_dekompozicija(P);
        A_cpy = A_cpy.primijeni_permutaciju(P);

        // Rješavanje sustava za svaki standardni bazni vektor za dobivanje stupaca inverzne matrice
        for (int i = 0; i < n; ++i) {
            Matrica e(n, 1);
            e[i][0] = 1;

            Matrica tmp = A_cpy.kopiraj();
            Matrica Pe = e.primijeni_permutaciju(P);

            // Supstitucije unaprijed i unatrag za svaki stupac
            tmp.supstitucija_unaprijed(Pe);
            tmp.supstitucija_unatrag(Pe);

            // Postavljanje rješenja kao stupca inverzne matrice
            A_inv.postavi_stupac(i, Pe);
        }
        return A_inv;
    }

    /**
     * Metoda izračunava determinantu matrice koristeći LUP dekompoziciju.
     * @return Determinanta matrice.
     */
    double det() {
        // Provjera je li matrica kvadratna
        if (rows != cols) {
            throw std::invalid_argument("Dimenzije matrice i vektora nisu kompatibilne za računanje determinante.");
        }

        double det;
        int swap;
        int n = rows;
        Matrica P(n, 1);

        // Izvršavanje LUP dekompozicije i brojqnje permutacija
        swap = LUP_dekompozicija(P);
        *this = primijeni_permutaciju(P);

        // Postavljanje predznaka determinante na temelju parnosti broja permutacija
        if (swap % 2 != 0 ) {
            det = -1.;
        }
        else {
            det = 1.;
        }

        // Množenje elemenata na dijagonali gornje trokutaste matrice
        // Ako je dekompozicija obavljena, umnožak po dijagonali donje trokutaste matrice iznosi 1
        for (int i = 0; i < n; ++i) {
            det *= elements[i][i];
        }
        return det;
    }

};

// Koristi se za zapis rješenja zadatka u datoteku txt.
void zapisRjesenjaUDatoteku(const std::string& filename, const std::string& sadrzaj) {
    std::ofstream file(filename);
    if (file.is_open()) {
        file << sadrzaj;
        file.close();
    } else {
        std::cerr << "Ne mogu otvoriti datoteku za pisanje: " << filename << std::endl;
    }
}

void prvi() {
    std::cout << "1. zadatak" << std::endl;
    Matrica A("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/1/A.txt");
    Matrica B = A * 128.0;
    Matrica C = B * (1.0 / 128.0);

    bool jednake = (A == C);
    std::string rezultat = "Jesu li matrice jednake unutar epsilon tolerancije? " + std::string(jednake ? "Da" : "Ne") + "\n";
    std::cout << rezultat;

    zapisRjesenjaUDatoteku("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/1/rjesenje.txt", rezultat);
}

void drugi() {
    std::cout << "2. zadatak" << std::endl;
    Matrica A_lu("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/2/A.txt");
    Matrica b_lu("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/2/b.txt");
    Matrica A_lup("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/2/A.txt");
    Matrica b_lup("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/2/b.txt");

    try {
        Matrica x = A_lu.rijesi_sustav_LU(b_lu);
        std::cout << "Rj. sustava - LU dekompozicija: " << std::endl;
        x.ispisi();
        x.ispisiUDatoteku("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/2/rjesenje_LU.txt");
    } catch (const std::runtime_error& e) {
        std::cout << "Sustav nije moguće riješiti LU dekompozicijom." << std::endl;
    }

    try {
        Matrica x1 = A_lup.rijesi_sustav_LUP(b_lup);
        std::cout << "Rj. sustava - LUP dekompozicija: " << std::endl;
        x1.ispisi();
        x1.ispisiUDatoteku("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/2/rjesenje_LUP.txt");
    } catch (const std::runtime_error& e) {
        std::cout << "Sustav nije moguće riješiti LUP dekompozicijom." << std::endl;
    }
}

void treci() {
    std::cout << "3. zadatak" << std::endl;
    Matrica A_lu("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/3/A.txt");
    Matrica A_lup("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/3/A.txt");

    Matrica A_lu_sustav("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/3/A.txt");
    Matrica b_lu("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/3/b.txt");
    Matrica A_lup_sustav("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/3/A.txt");
    Matrica b_lup("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/3/b.txt");

    int n = A_lu.rows;
    Matrica P(n, 1);
    A_lu.LU_dekompozicija();
    A_lup.LUP_dekompozicija(P);
    A_lup = A_lup.primijeni_permutaciju(P);

    std::cout << "Matrica za LU dekompoziciju: " << std::endl;
    A_lu.ispisi();
    std::cout << "Matrica za LUP dekompoziciju: " << std::endl;
    A_lup.ispisi();
    A_lu.ispisiUDatoteku("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/3/matrica_lu.txt");
    A_lup.ispisiUDatoteku("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/3/matrica_lup.txt");

    try {
        Matrica x = A_lu_sustav.rijesi_sustav_LU(b_lu);
        x.ispisi();
        x.ispisiUDatoteku("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/3/rjesenje_sustava_LU.txt");
    } catch (const std::runtime_error& e) {
        std::string rezultat = "Sustav nije moguće riješiti LU dekompozicijom.\n";
        std::cout << rezultat;
        zapisRjesenjaUDatoteku("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/3/rjesenje_sustava_LU.txt", rezultat);
    }
    try {
        Matrica x = A_lup_sustav.rijesi_sustav_LUP(b_lup);
        x.ispisi();
        x.ispisiUDatoteku("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/3/rjesenje_sustava_LUP.txt");
    } catch (const std::runtime_error& e) {
        std::string rezultat = "Sustav nije moguće riješiti LUP dekompozicijom.\n";
        std::cout << rezultat;
        zapisRjesenjaUDatoteku("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/3/rjesenje_sustava_LUP.txt", rezultat);
    }

}

void cetvrti() {
    std::cout << "4. zadatak" << std::endl;
    Matrica A_lu("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/4/A.txt");
    Matrica b_lu("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/4/b.txt");
    Matrica A_lup("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/4/A.txt");
    Matrica b_lup("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/4/b.txt");

    try {
        Matrica x = A_lu.rijesi_sustav_LU(b_lu);
        x.ispisi();
        x.ispisiUDatoteku("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/4/rjesenje_LU.txt");
    } catch (const std::runtime_error& e) {
        std::cout << "Sustav nije moguće riješiti." << std::endl;
    }

    try {
        Matrica x1 = A_lup.rijesi_sustav_LUP(b_lup);
        x1.ispisi();
        x1.ispisiUDatoteku("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/4/rjesenje_LUP.txt");
    } catch (const std::runtime_error& e) {
        std::cout << "Sustav nije moguće riješiti." << std::endl;
    }
}

void peti() {
    std::cout << "5. zadatak" << std::endl;
    Matrica A("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/5/A.txt");
    Matrica b("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/5/b.txt");

    try {
        Matrica x = A.rijesi_sustav_LUP(b);
        x.ispisi();
        x.ispisiUDatoteku("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/5/rjesenje.txt");
    } catch (const std::runtime_error& e) {
        std::cout << "Sustav nije moguće riješiti." << std::endl;
    }
}

void sesti() {
    std::cout << "6. zadatak" << std::endl;
    Matrica A("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/6/A.txt");
    Matrica b("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/6/b.txt");
    EPSILON = 10e-10;

    try {
        Matrica x = A.rijesi_sustav_LUP(b);
        x.ispisi();
        x.ispisiUDatoteku("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/6/rjesenje.txt");
    } catch (const std::runtime_error& e) {
        std::cout << "Sustav nije moguće riješiti LUP dekompozicijom." << std::endl;
    }

    EPSILON = 10e-8;
}

void sedmi() {
    std::cout << "7. zadatak" << std::endl;
    Matrica A("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/7/A.txt");

    try {
        Matrica A_inv = A.inverz();
        A_inv.ispisi();
        A_inv.ispisiUDatoteku("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/7/rjesenje.txt");
    } catch (const std::runtime_error& e) {
        std::string rezultat = "Inverz za matricu ne postoji, singularna matrica.\n";
        std::cout << rezultat;
        zapisRjesenjaUDatoteku("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/7/rjesenje.txt", rezultat);
    }
}

void osmi() {
    std::cout << "8. zadatak" << std::endl;
    Matrica A("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/8/A.txt");

    try {
        Matrica A_inv = A.inverz();
        A_inv.ispisi();
        A_inv.ispisiUDatoteku("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/8/rjesenje.txt");
    } catch (const std::runtime_error& e) {
        std::cout << "Inverz za matricu ne postoji, singularna matrica." << std::endl;
    }
}

void deveti() {
    std::cout << "9. zadatak" << std::endl;
    Matrica A("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/9/A.txt");

    try {
        double det = A.det();
        std::string rezultat = "Det(A) = " + std::to_string(det) + "\n";
        std::cout << rezultat;
        zapisRjesenjaUDatoteku("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/9/rjesenje.txt", rezultat);
    } catch (const std::invalid_argument& e) {
        std::cout << "Matrica nije kvadratna." << std::endl;
    }
}

void deseti() {
    std::cout << "10. zadatak" << std::endl;
    Matrica A("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/10/A.txt");

    try {
        double det = A.det();
        std::string rezultat = "Det(A) = " + std::to_string(det) + "\n";
        std::cout << rezultat;
        zapisRjesenjaUDatoteku("/Users/dominik/Desktop/Diplomski/1semestar/APR/0036538320/10/rjesenje.txt", rezultat);
    } catch (const std::invalid_argument& e) {
        std::cout << "Matrica nije kvadratna." << std::endl;
    }
}

int main() {
    prvi();
    drugi();
    treci();
    cetvrti();
    peti();
    sesti();
    sedmi();
    osmi();
    deveti();
    deseti();
    return 0;
}
