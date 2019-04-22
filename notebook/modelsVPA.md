Modèles pour Vlasov 1dx1dv
====

Le champ électrique, nécessaire dans les différentes modélisations, peut être traité soit à l'aide de l'équation de Poisson :

$$
  \partial_x E = \int f\,\mathrm{d}v - 1 = \rho - 1
$$

soit à l'aide de l'équation d'Ampère :

$$
  \partial_t E = -\int vf\,\mathrm{d}v = - j
$$

# Modélisation cinétique (K)

$$
  \begin{cases}
    \partial_t f + v\partial_x f + E\partial_v f = 0 \\
    \partial_x E = \int f\,\mathrm{d}v \quad \text{ou} \quad \partial_t E = -\int vf\,\mathrm{d}v
  \end{cases}
$$

# Modèle hybride

On cherche à représenter une distribution de particules chaudes et de particules froides. La dynamique des particules froides étant proche d'un état d'équilibre, elles seront modélisée à l'aide d'un modèle fluide, avec les équations d'Euler. On se place dans le cas extrême sans température, la distribution en vitesse de telles particules est représentée par une distribution de Dirac.

$$
  f = f_c + f_h \quad , \quad f_c = n_c(t,x)\delta_{u_c(t,x)}(v)
$$

## Modèle hybride non-linéarisé (HNL)

$$
  \begin{cases}
    \partial_t n_c + \partial_x(n_cu_c) = 0 \\
    \partial_t(n_cu_c) + \partial_x(n_cu_c^2) = n_cE \\
    \partial_t f_h + v\partial_x f_h + E\partial_v f_h = 0 \\
    \left( \partial_x E = n_c + \int f_h\,\mathrm{d}v - 1 \right) \\
    \partial_t E = -n_cu_c - \int vf_h\,\mathrm{d}v
  \end{cases}
$$

## Modèle hybride linéarisé (HL)

On souhaite linéarisé autour de l'état :

$$
  n_c = n_c^0(t,x) + \epsilon\tilde{n}_c \quad u_c = \epsilon\tilde{u}_c \quad E = E^0 + \epsilon\tilde{E}
$$

La première équation d'Euler nous donne :

$$
  \partial_t(n_c^0 + \epsilon\tilde{n}_c) + \epsilon\partial_x(n_c^0\tilde{u}_c) = 0 \Rightarrow n_c^0 \perp\!\!\!\perp t
$$

Ce qui permet de simplifier cette équation :

$$
  \partial_t\tilde{n}_c + \partial_x(n_c^0\tilde{u}_c) = 0
$$

L'équation de Poisson nous donne :

$$
  \partial_x (E_0+\epsilon\tilde{E}) = n_c^0 + \epsilon\tilde{n}_c + \int f_h^0\,\mathrm{d}v + \epsilon\int \tilde{f}_h\,\mathrm{d}v - n_e
$$

En ne considérant que les termes d'ordre 1 en $\epsilon$, on trouve :

$$
  \partial_x E_0 = n_c^0 + \int f_h^0\,\mathrm{d}v - n_e
$$

On peut imposer à la condition initiale qui correspond à un équilibre $n_c^0 + \int f_h^0\,\mathrm{d}v - n_e = 0$, notre plasma est quasi-neutre, soit $E_0 = 0$. Il s'agit là d'une hypothèse physique.

On obtient alors :

$$
  \begin{cases}
    \partial_t \tilde{n}_c + \partial_x(n_c^0\tilde{u}_c) = 0 \\
    \partial_t \tilde{u}_c = \tilde{E} \\
    \partial_t f_h + v\partial_xf_h + \epsilon\tilde{E}\partial_vf_h = 0 \\
    \epsilon \partial_t \tilde{E} = -\epsilon n_c^0\tilde{u}_c - \int vf_h\,\mathrm{d}v
  \end{cases}
$$

La première équation est la seule faisant intervenir $\tilde{n}_c$, elle n'est donc pas nécessaire pour obtenir les autres variables, on ne s'intéresse donc qu'au modèle réduit :

$$
  \begin{cases}
    \partial_t \tilde{u}_c = \tilde{E} \\
    \partial_t f_h + v\partial_xf_h + \epsilon\tilde{E}\partial_vf_h = 0 \\
    \epsilon \partial_t \tilde{E} = -\epsilon n_c^0\tilde{u}_c - \int vf_h\,\mathrm{d}v
  \end{cases}
$$ 

Grâce à l'hypothèse de quasi-neutralité initiale, on peut réécrire l'équation de Poisson :

$$
  \epsilon\partial_x\tilde{E} = n_c^0 + \epsilon\tilde{n}_c + n_h^0 + \epsilon\tilde{n}_h - n_e
$$

soit :

$$
  \partial_x\tilde{E} = \tilde{n}_c + \tilde{n}_h
$$

avec $\tilde{n}_h = \int (f_h-f_h^0)\,\mathrm{d}v$

L'équation d'Ampère nous donne :

$$
  \epsilon\partial_t\tilde{E} = -\epsilon n_c^0 \tilde{u}_c - \int v(f_h^0 + \epsilon\tilde{f}_h)\,\mathrm{d}v
$$

En ne considérant que les termes d'ordre 1 on obtient : $\int vf_h^0\,\mathrm{d}v = 0$, on peut alors écrire l'équation d'Ampère comme :

$$
  \partial_t\tilde{E} = n_c^0\tilde{u}_c - \tilde{\jmath}_h
$$

avec $\tilde{\jmath}_h = \int v(f_h - f_h^0)\,\mathrm{d}v$

Le modèle linéarisé, avec pour inconnue $\tilde{u}_c$, $\tilde{E}$ et $f_h$ nous donne :

$$
  \begin{cases}
    \partial_t\tilde{E} = -n_c^0\tilde{u}_c - \tilde{\jmath}_c \\
    \partial_x\tilde{E} = \tilde{n}_c \tilde{n}_h \\
    \partial_t\tilde{u}_c = \tilde{E} \\
    \partial_tf_h + v\partial_xf_h + \epsilon\tilde{E}\partial_vf_h = 0
  \end{cases}
$$

avec :

$$
  \tilde{\jmath}_h = \int v(f_h - f_h^0)\,\mathrm{d}v \qquad \tilde{n}_h = \int (f_h - f_h^0)\,\mathrm{d}v
$$

# Conditions initiales

Pour effectuer une simulation et comparer les modèles, il est nécessaire d'avoir des conditions initiales similaires pour les modèles. Les modèles sont ceux qui se prêtent le mieux à l'écriture de la condition initiale comme une bi-maxwellienne symétrique (particules chaudes) autour d'une distribution de Dirac (des particules froides).

Pour le modèle hybride (linéarisé ou non) la condition initiale est :

* $$f_h^0 = \frac{n_h^0}{2\sqrt{2\pi}}\left[ \exp{\left(-\frac{|v-u_h^0|^2}{2}\right)} + \exp{\left(-\frac{|v+u_h^0|^2}{2}\right)} \right]$$
* $n_c^0 = 1$ et on doit vérifier $n_h^0 \ll n_c^0$
* $\tilde{u}_c(t=0) \approx 0 \approx \tilde{n}_c(t=0)$
* $f_h(t=0,x,v) = f_h^0(v)(1+\varepsilon\cos(kx))$

Pour le modèle cinétique, il n'est pas possible numériquement d'avoir une distribution de Dirac :

* $$f(t=0,x,v) = \frac{n_c^0}{\sqrt{2\pi}}\exp\left(-\frac{|v|^2}{2T_c}\right) + f_h(t=0)$$
* En vérifiant bien $T_c\ll 1$

Modèles pour Vlasov 1dx2dv
====

# Modélisation cinétique (K)

$$
  \begin{cases}
    \partial_t f + v\partial_x f + (E + v\times B)\partial_v f = 0 \\
    \partial_x E = \int f\,\mathrm{d}v \quad \text{ou} \quad \partial_t E = -\int vf\,\mathrm{d}v
  \end{cases}
$$



