---
title: Contre-exemples en simple précision
subtitle: Six défauts d’inclusion, dont un retard négatif sans fast-math
author: Yann Orlarey
date: 8 octobre 2026
document-style: article-a4
language: fr
---

::: toc+
- **Cadre des expériences** — définir le format, les domaines et le sens des assertions.
- **Un retard négatif sans fast-math** — montrer qu’une analyse non négative peut exclure un retard égal à −1.
- **Les autres contre-exemples** — isoler les défauts de conversion, d’arithmétique, de contraction et de réduction trigonométrique.
  - **Addition au-delà de 2²⁴** — comparer le point annoncé avec l’arrondi binary32.
  - **Multiplication sur un domaine non ponctuel** — observer la valeur située sous la borne inférieure.
  - **Conversion explicite d’un entier en float** — vérifier le résultat de FloatCast.
  - **Contraction FMA sans réassociation** — distinguer deux arrondis d’un arrondi fusionné.
  - **Sinus d’un argument modeste** — montrer une erreur de réduction d’argument, avec compensation activée.
  - **Diagnostic supplémentaire sur un grand argument** — isoler un comportement indéfini du calcul de précision.
- **Exécuter les tests** — reproduire les échecs d’inclusion et les vérifier avec l’oracle facultatif.
- **Ce que ces tests établissent** — préciser la portée des preuves et leur utilisation après correction.
:::

# Cadre des expériences

Les témoins utilisent directement la bibliothèque de ce dépôt, avec `itv::programPrecision() = 1` et `itv::libmCompensation() = true`. Ils ont été construits à partir du commit `7db5ec3`, sur la branche `fix/interval-bounds-safety`. Ils n’utilisent pas un binaire du compilateur Faust : l’acceptation d’un DSP particulier par Faust n’est pas supposée.

Le contrat du calcul témoin est binary32, avec arrondi au plus proche et départage des égalités vers le significand pair. Le code est compilé sans fast-math et avec contraction implicite désactivée. Des variables `volatile float` matérialisent chaque étape. Le témoin FMA utilise explicitement `std::fma` pour tester un contrat permettant cette contraction.

Les entrées des opérations flottantes sont des domaines dont les bornes sont déjà exactement représentables en float, construits avec un LSB négatif. Chaque valeur exécutée appartient à son domaine d’entrée. Le défaut ne vient donc ni d’une entrée hors contrat, ni d’un choix volontaire de l’arithmétique entière. Le cas `FloatCast` part, lui, d’un entier défini dans le domaine donné à la conversion.

Les six assertions demandent la même propriété :

```math
v_{\mathrm{exécuté}} \in I_{\mathrm{annoncé}}
```

**Les six assertions échouent actuellement.** Les résultats ont été reproduits avec Apple Clang 21, GCC 15, les détecteurs de comportements indéfinis de Clang et Emscripten 3.1.74 sous Node. L’oracle MPFR confirme séparément les valeurs binary32 de ces témoins.

| Test | Intervalle annoncé | Valeur exécutée |
|---|---|---:|
| `addition` | `[16777217, 16777217]` | `16777216` |
| `multiplication` | `[16785409, 16793604]` | `16785408` |
| `delay` | `[0, 8195]` | `−1` |
| `conversion` | `[16777217, 16777217]` | `16777216` |
| `contraction` | `[0, 0]` | `−1` |
| `sine` | `[-0.50636541843414307, -0.50636541843414307]` | `-0.50636565685272217` |

# Un retard négatif sans fast-math

Considérons les domaines de deux paramètres flottants indépendants :

```math
A = B = [4097, 4098]
```

L’expression analysée, avec cet ordre d’opérations, est :

```cpp
int((a * b - 16785408.0f) - 1.0f)
```

Le calcul réel exact est non négatif sur tout le domaine :

```math
4097^2 - 16785408 - 1 = 16785409 - 16785408 - 1 = 0
```

La bibliothèque annonce effectivement `[0, 8195]`. Le calcul float séparé produit pourtant `−1` à la borne inférieure :

| Étape | Analyse actuelle | Calcul float pour `a = b = 4097` |
|---|---|---:|
| `a * b` | `[16785409, 16793604]` | `16785408` |
| `produit - 16785408` | `[1, 8196]` | `0` |
| `différence - 1` | `[0, 8195]` | `−1` |
| conversion en entier | `[0, 8195]` | `−1` |

À cette magnitude, les floats sont espacés de `2`. Le produit exact `16785409` est à mi-chemin entre `16785408` et `16785410`. Le départage vers le significand pair choisit `16785408`.

La borne calculée par la bibliothèque reste `16785409`. Dans [programBound](interval/interval_def.hh), une exception conserve certaines grandes bornes de valeur entière au lieu de les convertir en float. Cette exception s’applique ici à un **résultat flottant** dont la valeur réelle exacte est entière. Elle conserve une borne non représentable en binary32 et exclut le produit exécuté.

Le témoin scalaire impose les trois arrondis :

```cpp
volatile float a = 4097.0f;
volatile float b = 4097.0f;
volatile float product = a * b;
volatile float difference = product - 16785408.0f;
volatile float value = difference - 1.0f;
int delay = int(value); // value = -1, conversion définie
```

Toutes les constantes sont exactement représentables en float. Il n’y a ni libm, ni FMA, ni réassociation, ni sous-normal, ni dépassement, ni conversion indéfinie. Le domaine de départ n’est pas réduit à un point.

L’intervalle des résultats de cette évaluation binary32 séparée est `[-1, 8195]`, par monotonie des trois opérations sur ces domaines positifs. Annoncer `[0, 8195]` exclut donc une extrémité effectivement atteinte. Une décision fondée sur la non-négativité du retard serait injustifiée.

# Les autres contre-exemples

## Addition au-delà de 2²⁴

Pour deux paramètres de valeurs fixées `a = 16777216.0f` et `b = 1.0f`, l’analyse annonce le point `16777217`. Le calcul float renvoie `16777216` : au-dessus de `2²⁴`, l’espacement vaut `2`, et l’égalité est départagée vers le significand pair.

La borne annoncée n’est pas un float. Le défaut apparaît dès une addition, avant toute composition ou conversion entière. Il provient de la même exception de [programBound](interval/interval_def.hh).

## Multiplication sur un domaine non ponctuel

Le test `multiplication` isole la première étape du retard précédent. Il demande que le résultat de `4097.0f * 4097.0f` appartienne à l’image analysée de `[4097, 4098] × [4097, 4098]`.

La valeur `16785408` se situe sous la borne `16785409`. Ce test restera pertinent si le traitement des constantes évolue : les deux facteurs sont des paramètres dans des domaines non ponctuels.

## Conversion explicite d’un entier en float

Le test `conversion` donne à `FloatCast` le résultat de `IntNum(16777217)`. Cet entier est exactement représentable en int32 et en double. Sa conversion binary32 vaut `16777216`.

`FloatCast` annonce pourtant `[16777217, 16777217]`. Forcer un LSB négatif ne suffit pas : l’exception sur les grandes bornes entières conserve encore la valeur non convertible exactement en float. Ce test concerne une conversion explicite, sans ambiguïté sur le type de l’opération.

## Contraction FMA sans réassociation

Les paramètres valent :

```math
\begin{aligned}
a &= 1 + 2^{-13} \\
b &= 1 - 2^{-13} \\
ab &= 1 - 2^{-26}
\end{aligned}
```

L’évaluation séparée arrondit le produit à `1`, puis obtient `0` après soustraction de `1`. La bibliothèque retrouve le point `[0, 0]` pour `a*b - 1`.

Une FMA fusionnant uniquement le produit et cette soustraction conserve leur résultat exact `−2⁻²⁶`. En multipliant ensuite par `2²⁶`, on obtient un entier égal à `−1`, tandis que l’analyse conserve `[0, 0]`.

```cpp
volatile float a = 1.0f + 0x1p-13f;
volatile float b = 1.0f - 0x1p-13f;
volatile float fused = std::fma(float(a), float(b), -1.0f);
volatile float value = fused * 0x1p26f;
int delay = int(value); // -1
```

**Ce cas suppose un contrat autorisant la contraction.** Il ne prétend pas qu’une exécution interdisant toute FMA produirait ce résultat. Aucune réassociation n’est nécessaire : la fusion concerne le produit et sa soustraction immédiate. La FMA est explicite dans le test, ce qui rend la reproduction indépendante des options par défaut du compilateur.

Un futur encadrement conservateur pour le contrat « FMA permise, sans réassociation » devra couvrir cette valeur. Une analyse reproduisant uniquement les arrondis séparés ne le fait pas.

## Sinus d’un argument modeste

L’entrée est le float exactement représentable `100.0f`. La bibliothèque annonce `[-0.50636541843414307, -0.50636541843414307]`. Le sinus float exécuté vaut `-0.50636565685272217`, et l’oracle MPFR, arrondi en binary32, confirme cette dernière valeur dans les environnements testés.

Le problème vient de [SinBounds](interval/intervalSin.cpp) : le mode float réduit l’argument avec `fmod(x, 2*M_PI)`, où `M_PI` est une approximation double de π, puis construit un intervalle qui arrondit la phase réduite en float. Le sinus de cette phase arrondie ne reproduit pas le sinus de l’entrée initiale. Une valeur ponctuelle obtenue dans cet intermédiaire n’est pas une preuve que cette réduction conserve la valeur recherchée.

La compensation libm est activée dans le test. Le résultat ponctuel est actuellement exempt de cette compensation en mode float. Les deux valeurs diffèrent de quatre ULPs binary32 à cette magnitude ; une marge de deux ULPs ne suffirait pas ici. L’erreur concerne aussi le calcul de référence de l’analyseur, pas seulement l’arrondi de la libm cible.

La suite compare la valeur de la libm réellement exécutée. Sa variante MPFR compare en plus le sinus mathématique correctement arrondi en binary32, sans présumer que toute libm fournit universellement ce même arrondi.

## Diagnostic supplémentaire sur un grand argument

Le cas optionnel `sine-large` utilise `1.0e20f`, dont la valeur exacte est `100000002004087734272`. Il ne fait pas partie des six assertions précédentes, car il révèle **un comportement indéfini dans l’analyseur**.

Dans [exactPrecisionUnary](interval/precision_utils.hh), l’évaluation de `sin(x + u) - sin(x)` donne zéro pour ce grand argument et le petit pas utilisé. Le calcul suivant produit `log2(0) = -∞`, puis tente une conversion en `int`, qui n’est pas définie. La vérification ultérieure de `INT_MIN` dans `SinBounds` arrive trop tard pour rendre cette conversion valide.

Dans une construction instrumentée, le diagnostic s’arrête avec :

```text
precision_utils.hh:62:15: runtime error:
-inf is outside the range of representable values of type 'int'
```

Sans instrumentation, les constructions natives testées affichent un intervalle ponctuel proche de `−0.17085628`, contre un sinus float et une référence MPFR proches de `0.65657669`. Ce résultat illustre en outre la faiblesse de la réduction avec un `2π` approché pour les grands arguments. **Il ne constitue pas une preuve d’exécution C++ définie**, contrairement aux six témoins principaux : l’analyseur a déjà rencontré un comportement indéfini.

Le diagnostic reste exécutable séparément et n’est pas masqué dans les tests :

```sh
./build-float-witnesses/FloatConservativenessTests --case sine-large
```

Pour reproduire le diagnostic de conversion :

```sh
cmake -S . -B build-float-ubsan -DNOTIDY=ON \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-fsanitize=undefined,float-cast-overflow -fno-sanitize-recover=all"
cmake --build build-float-ubsan --parallel
./build-float-ubsan/FloatConservativenessTests --case sine-large
```

Les six cas principaux passent leurs contrôles de validité avec cette instrumentation ; leurs assertions d’inclusion échouent toujours. Le cas supplémentaire s’arrête sur le comportement indéfini.

# Exécuter les tests

Les sources sont dans [float_conservativeness_tests.cpp](tests/float_conservativeness_tests.cpp). Construire sans dépendance numérique extérieure :

```sh
cmake -S . -B build-float-witnesses -DNOTIDY=ON
cmake --build build-float-witnesses --parallel
./build-float-witnesses/FloatConservativenessTests --require-inclusion
```

La dernière commande retourne actuellement **le code `1` et six échecs d’inclusion**. Chaque échec imprime l’intervalle annoncé et la valeur exclue. Un cas peut être exécuté seul :

```sh
./build-float-witnesses/FloatConservativenessTests --case delay
```

La même source fournit un oracle facultatif :

```sh
cmake -S . -B build-float-oracle -DNOTIDY=ON \
  -DINTERVAL_ENABLE_MPFR_TESTS=ON
cmake --build build-float-oracle --parallel
./build-float-oracle/FloatConservativenessOracleTests --require-inclusion
```

Cette dernière commande retourne également `1`. MPFR utilise une précision de 24 bits et arrondit **après chaque opération du programme témoin**. Tous les résultats intermédiaires des six cas sont finis, normaux ou nuls ; aucun ne nécessite un arrondi sous-normal. Il n’est donc pas nécessaire d’émuler ici le sous-dépassement binary32. Ce dispositif n’est pas encore un oracle général couvrant toutes les particularités du float.

Pour WebAssembly, après activation d’Emscripten :

```sh
emcmake cmake -S . -B build-wasm-float -DNOTIDY=ON \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build-wasm-float --parallel
node build-wasm-float/FloatConservativenessTests.js --require-inclusion
```

CTest emploie un mode explicitement différent :

```sh
ctest --test-dir build-float-witnesses -L known-float-gaps -V
```

Les tests enregistrés passent `--expect-known-gaps`. Leur succès signifie **que les défauts connus ont été reproduits**, et non que l’inclusion float est garantie. Une correction qui couvre une valeur fait échouer ce mode : il faudra alors transférer le témoin corrigé dans la suite exigeant l’inclusion. Un témoin invalide ou une erreur de paramètres reste une erreur, même dans le mode des défauts connus.

# Ce que ces tests établissent

L’absence de conservation est démontrée par des valeurs numériques exclues, avec des entrées respectant leurs domaines. Le témoin de retard montre une conséquence après conversion entière, dans une exécution stricte sans libm ni FMA.

Les causes sont distinctes : conservation d’une borne entière non représentable alors que l’opération est flottante, contraction autorisée non couverte par les arrondis séparés, et réduction trigonométrique dont l’arrondi intermédiaire n’est pas encadré. Le diagnostic optionnel ajoute un défaut de conversion dans les métadonnées de précision. Les tests de régression float précédents pouvaient passer tout en laissant ces défauts présents.

Les corrections devront faire passer les assertions d’inclusion ordinaires. Le passage de ces six témoins sera une condition nécessaire ; une garantie générale exigera encore des règles conservatrices pour toutes les opérations et un contrat couvrant les exécutions analysées.
