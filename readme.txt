# Upute za pokretanje programa

Ovaj projekt je razvijen kao dio domaće zadaće za predmet APR.
Program sadrži funkcije za rad s matricama, uključujući operacije poput množenja, zbrajanja, oduzimanja, dekompozicija i rješavanje sustava linearnih jednadžbi.

    ## Preduvjeti

    Za uspješno prevođenje i pokretanje programa potrebno je sljedeće:
    - Razvojno okruženje koje podržava CMake (npr. CLion, Visual Studio Code)
    - Instaliran C++ prevoditelj (npr. GCC ili Clang)
    - CMake alat (verzija 3.29 ili novija)

    ## Struktura projekta

    0036538320
    - CMakeLists.txt
    - main.cpp
    - 1/...
    - 2/...
    ...
    - 10/...
    - readme.txt

    ## Upute za prevođenje

        ### Korištenje CLiona
        1. Otvori CLion i odaberi `File -> Open`, zatim pronađi i odaberi direktorij `0036538320`.
        2. Kada se projekt učita, odaberi `Build -> Build Project` za prevođenje cijelog projekta.
        3. Ako postoje greške pri prevođenju, provjeri putanje ulaznih datoteka i postavke projekta.

        ### Prevođenje iz terminala
        Ako koristiš terminal, možeš prevesti projekt koristeći sljedeće naredbe:
        1. Otvori terminal i idi do direktorija `0036538320`.
        2. Pokreni sljedeće naredbe:
            ```
            mkdir build
            cd build
            cmake ..
            make
            ```
        3. Ovo će stvoriti izvršnu datoteku u `build` direktoriju.

    ## Upute za pokretanje programa

    Nakon što je projekt uspješno preveden, možeš pokrenuti program na sljedeći način:

        ### Pokretanje u CLionu
        - Koristi `Run` opciju u CLionu za pokretanje programa.
        - Program će izvršiti zadatke i ispisati rezultate na konzoli te u odgovarajuće `.txt` datoteke.

        ### Pokretanje iz terminala
        1. Prevedeni izvršni program nalazi se u `build` direktoriju. Pokreni ga s:
            ```
            ./ime_izvršne_datoteke
            ```
        2. Program će pročitati ulazne datoteke, izvršiti odgovarajuće zadatke i zapisati rezultate u izlazne `.txt` datoteke.

    ## Dodatne napomene

    - Ulazne datoteke (`A.txt`, `b.txt` itd.) nalaze se u odgovarajućim poddirektorijima (`1`, `2`...) unutar direktorija projekta.
    - Rezultati izvođenja zadataka bit će zapisani u odgovarajuće izlazne `.txt` datoteke unutar tih istih direktorija.
