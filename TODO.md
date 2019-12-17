# TODO


Liste des choses à faire prochainement :

- [x] Transport par WENO en C++ (utiliser la classe `multi_array` de `boost`)
- [x] Transport par FFT en C++ (utiliser `fftw`)
- [x] Écrire le schéma IFRKSSP(3,3) et IFRKSSP(4,3) correspondant
    > En phase de réalisation (schéma IFRK(4,4) plus ou moins écrit), avec un script python pour s'assurer de la validité des résultats

- [ ] Passer à VP 1dx2dv
- [ ] Trouver un nom pour la simulation par FFT+WENO avec un schéma Lawson
    + ~~**F**FT, **W**ENO, **V**lasov **P**oisson : `f3vp`~~
    + ~~F**F**T, W**E**NO, L**a**wson : `fea` (il semblerait qu'en espagnol, cela signifie *moche*, donc pas terrible)~~
    + **S**pectral W**E**NO La**w**son : `sew` (coudre c'est sympa comme nom, ça permet de *filer* la métaphore en parlant de fil, de broderie ou je ne sais quoi).
    + WENO FFT/Spectral IFRK/Lawson Vlasov Poisson
- [ ] Implémenter un solveur d'Euler 1D avec WENO en C++
    + [ ] Vérifier le code en régime fluide
    + [ ] Tester sur du hybride fluide-cinétique

Il est aussi envisageable d'effectuer des tests unitaires :

- [ ] Méthdes du WENO (tester sur une itération la valeur de l'erreur inférieur à une certaine valeur, vérifier l'ordre, ou autre)
- [ ] Transport par FFT (compliqué de tester la résolution de l'équation de Poisson, mais possibilité de tester un transport)
- [ ] Tester le calcul de densité et autres opérations sur un `field<_T>` sur des cas informatiquement simple (somme des n premiers entiers, etc.)

Échéancier en temps court :

- [x] Tester avec des `std::complex<double>` et comparer les perfs (à mon avis effet non significatif sur WENO+FFT, et Strang est déjà assez rapide)
    > La fonction `std::exp(std::complex<_T> const&)` renvoie un `std::polar<_T>` ce qui lui donne un coût très réduit. Il y a ensuite la conversion de `std::polar<_T>` en `std::complex<_T>` puis le cast vers `fftw_complex`. Le gain en facilité d'écriture est tel que je pense que cela vaut le coup/coût.
    
    > Fait dans la classe `spectrum_` (le temps de convertir tout le code avec ça). Pas de comparaison de perfs effectuées mais la simplicité d'écriture l'emporte largement. (surtout pour WENO+FFT)
     
- [x] Tester avec IFRK(4,4), et potentiellement d'autres (écrire un script python pour l'écriture automatique d'un schéma IFRK à partir d'un schéma RK)
    - [x] Faire un script python qui génère le schéma IFRK associé à un schéma RK (par tableau de Butcher)
    - [x] Simplifier le schéma ainsi obtenu pour limiter le nombre d'appels à la fonction $N$

- [x] Tracer l'énergie au cours du temps : $$H(t) = \int v^2f\,\mathrm{d}v\mathrm{d}x + \int E^2\,\mathrm{d}x$$ Normalement c'est une constante. Tester ceci avec plusieurs schémas en temps.
    > Fait avec la fonction `field.h: _T energy(field<_T,NumDimsV> const&,ublas::vector<_T> const& E)` (ce qui me plaît pas est que j'impose un conteneur pour l'énergie $E$, sinon ça fonctionne. Testé que sur le code `strang.cc` où l'erreur relative évolue de manière similaire à l'énergie électrique, avec une différence de 4%.
    
    **TODO:**
    
    - [x] Savoir comment mesurer $(H(t^n))_n$, puis tracer $\mathcal{E}(\Delta t) = \left(\frac{H(t^n)-H(0)}{H(0)}\right)_{\Delta t} = \mu((H(t^n))_n)$
        > Utilisation de la norme $L^{\infty}$
        
    - [x] Tester avec le code WENO+FFT
        > Résultats compatibles avec l'énergie électrique, et la diffusion du schéma WENO
        
    - [x] Comparer la courbe $\mathcal{E}(\Delta t)$ pour différents schémas de Lawson
        > Fait au moins pour IFRK(3,3) et IFRK(4,4)
    
- [x] Valider l'ordre de la méthode de Strang par rotation ou translation d'un truc dont on connait la solution exacte.
    > Fait avec une rotation : ordre 2
    
- [x] Demander à Nicolas à propos du bruit de la méthode de Strang
    > Résultats normaux, compatible avec une FFT. Valeur intégrale non bruitée

