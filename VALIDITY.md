---
title: Validité des calculs d’intervalles
subtitle: Introduction de mayBeInvalid et obligations d’intégration dans Faust
date: 8 octobre 2026
document-style: article-a4
language: fr
---

::: toc+
- **Contrat et représentation** — distinguer les résultats numériques et les exécutions potentiellement invalides.
- **Règles de propagation** — couvrir les domaines, les conversions et les récurrences sans perdre les alertes.
- **Tests** — vérifier les deux précisions et les deux couches avec des témoins indépendants.
- **Intégration dans Faust** — préserver le flag et exiger un diagnostic lorsque la sûreté n’est pas établie.
:::

# Contrat et représentation

Le [critère de correction normatif du README](README.md#correctness-criterion-normative)
reste la référence : l’intervalle implémenté doit contenir l’encadrement idéal
`Y`, et son flag ne doit jamais valoir `false` lorsque le flag idéal vaut `true`.
Un élargissement ou un signalement prudent d’invalidité est autorisé.

Les deux couches possèdent désormais cet attribut, indépendant de la nature
entière ou flottante :

| Couche | Accès au flag |
| :--- | :--- |
| Intervalle ordinaire | `x.mayBeInvalid()` |
| Intervalle affine | `x.mayBeInvalid` |

La convention `lsb >= 0 ⇔ entier` reste inchangée. Le flag signale un NaN
possible ou une opération potentiellement indéfinie, notamment une conversion
flottante hors plage, un reste entier interdit ou un décalage hors contrat.

Les états suivants sont distincts :

| Bornes numériques | Flag | Signification |
| :--- | :--- | :--- |
| Non vides | `false` | Résultats numériques sans invalidité possible dans le contrat |
| Non vides | `true` | Résultats numériques et invalidité possible |
| Vides | `false` | Aucun résultat ni exécution invalide représenté : bottom |
| Vides | `true` | Exécution potentiellement invalide sans résultat numérique représenté |

`empty()` et `aempty()` construisent bottom. Une injection de NaN crée un état
uniquement invalide. Les NaN stockés comme sentinelles des bornes vides ne sont
pas des valeurs numériques de l’intervalle.

Le constructeur ordinaire accepte un quatrième argument :

```cpp
itv::interval x(-1, 4, -24, false);
itv::interval_algebra algebra;
auto y = algebra.Sqrt(x);  // [0, 2], mayBeInvalid() == true
```

`withInvalid(true)` ajoute une alerte en conservant les bornes et le `lsb`.
`withInvalid(false)` ne peut pas effacer une alerte précédente. `isEmpty()`
porte uniquement sur les bornes numériques ; `isValid()` exige aussi un flag
à `false`. Les prédicats de singleton ordinaire ne permettent pas de plier un
point qui porte une alerte. Le prédicat affine `isConst()` décrit, lui, un
corridor horizontal ; il ne prouve ni un singleton ni l’absence d’invalidité.

# Règles de propagation

Les transferts publics ordinaires ajoutent la validité après le calcul de leurs
bornes. Leurs noyaux numériques sont internes. Cette séparation empêche un
retour anticipé, notamment vers bottom, de perdre un flag d’entrée. Les règles
publiques affines suivent le même principe et leurs opérations non linéaires
consultent l’oracle ordinaire.

Les contrôles de domaine portent sur les opérandes du programme après les
conversions implicites. Ils distinguent les opérations entières avec wrapping,
qui respectent le contrat, des combinaisons qui produisent un NaN ou exécutent
une opération indéfinie. Les domaines partiellement valides conservent leur
image numérique et portent une alerte.

Le contrat numérique courant admet les infinis IEEE lorsqu’ils ne sont pas
NaN : `log(0)` ou `1/0` produisent des valeurs numériques infinies, tandis que
`0/0`, `0 × ∞` et `∞ − ∞` sont invalides. Exiger des résultats audio finis
nécessite aussi un contrôle de finitude. Les zéros signés restent fusionnés dans
la représentation ; les réciproques et puissances négatives en float et double
incluent les deux signes d’infini lorsque cela est nécessaire.

Pour `IntCast`, **le résultat tronqué vers zéro doit être représentable**.
En double, `2147483647.75` donne légalement `2147483647`, mais `2147483648`
est exclu. En float, le premier nombre peut être arrondi à `2147483648` avant
le cast : c’est la valeur effectivement reçue qui décide de la validité.

Une entrée uniquement hors plage donne des bornes numériques vides avec une
alerte. Une entrée partiellement admissible conserve ses résultats tronqués
valides. Les opérations bitwise et les décalages ne perdent pas l’alerte d’une
conversion préalable. Le flag ne rend jamais défini un cast C++ interdit.

Le chemin des puissances entières avec wrapping exige un exposant non négatif.
Un exposant entier potentiellement négatif porte donc une alerte ; il ne doit
pas être interprété silencieusement comme zéro. Les puissances négatives
flottantes suivent leur propre domaine.

La réunion combine les flags par OU, y compris lorsqu’un opérande a des bornes
vides. Le raffinement numérique par intersection conserve les alertes
antérieures. L’égalité, l’inclusion et l’ordre affine observent le flag : des
bornes inchangées avec une nouvelle alerte représentent un nouvel état. Le
widening affine conserve les alertes anciennes et nouvelles. Les deux ponts
`fromItv` et `toItv` préservent le flag, y compris pour les états uniquement
invalides et les corridors entiers normalisés après wrapping.

Les sélecteurs, contrôles, effets, paramètres d’interface et échantillons d’une
waveform propagent leurs alertes. Une assertion de bornes ne certifie pas
silencieusement des valeurs exclues. Les fonctions étrangères sans contrat
restent non certifiées ; un stockage flottant étranger inconnu peut contenir
un NaN.

# Tests

La cible C++17 `InvalidityTests` couvre notamment :

- Les quatre états de représentation et la distinction entre NaN injecté et
  sentinelle de bottom.
- La propagation ordinaire et affine, pour chaque opérande des opérations
  unaires et binaires, ainsi que les retours sur des bornes vides.
- Les domaines complets, partiellement valides et uniquement invalides ; les
  arguments infinis et les combinaisons binaires interdites.
- Les conversions aux frontières int32, les valeurs fractionnaires encore
  admissibles après troncature, les valeurs float voisines et les conversions
  successives entier → float → entier.
- Les restes entiers interdits, les comptes de décalage et le wrapping valide.
- Les changements de flag sans changement de bornes, réunions, raffinements,
  widening et ponts affines.
- Les opérandes des contrôles et des interfaces, y compris les paramètres par
  défaut, ainsi que les incertitudes des ressources externes.
- Les deux indices de table du corpus Faust : le compteur borné dans la boucle
  et celui qui wrappe avant un modulo extérieur.

Les témoins d’exécution vérifient que chaque NaN observé exige une alerte.
Les tests ne supposent pas que les valeurs finies de la libm sont correctement
arrondies. Des motifs de bits flottants à graine fixe servent de témoins
indépendants pour les conversions : **aucun cast C++ invalide n’est exécuté
par le test**.

La cible facultative `InvalidityOracleTests` ajoute un oracle MPFR à 256 bits
pour les domaines unaires/binaires et la troncature vers int32. MPFR/GMP ne sont
pas des dépendances de production.

```sh
cmake -S . -B build -DNOTIDY=ON -DINTERVAL_ENABLE_MPFR_TESTS=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

Les suites complètes sont vérifiées avec Clang et GCC, avec les oracles
facultatifs, ainsi qu’avec AddressSanitizer/UndefinedBehaviorSanitizer et en
WebAssembly sous Node. Les suites existantes de bornes float/double restent
actives ; le changement de validité ne dispense pas de l’inclusion numérique.
Les tests complètent les justifications des règles et ne constituent pas une
preuve pour toutes les exécutions ou toutes les libm.

# Intégration dans Faust

**L’ajout dans cette bibliothèque ne modifie pas automatiquement le domaine
de types du compilateur ni sa politique de refus.** L’intégration doit :

1. Conserver le flag dans les constructions et copies de types, caches,
   sérialisations, réunions et récurrences. Copier seulement `lo`, `hi` et `lsb`
   perd l’information. Pour reconstruire des bornes vides, utiliser
   `empty(lsb).withInvalid(flag)` ; des bornes NaN injectées décrivent une valeur
   invalide et ne reconstruisent pas bottom.
2. Vérifier les domaines des opérations avant leur pliage constant. Une
   conversion invalide ne doit pas disparaître derrière un entier déjà plié.
3. Exiger un flag à `false`, une conversion définie et des bornes dans la plage
   admissible avant de certifier un retard ou un indice de table. Sinon, refuser
   avec un diagnostic, en distinguant le risque certain de l’incertitude.
4. Fournir les contrôles de domaine propres aux ressources : taille des tables,
   nombre de canaux et parties des soundfiles. Les attributs de valeurs ne
   conservent pas ces dimensions ; les lectures génériques signalent donc
   l’impossibilité de certifier. Les écritures ordinaires et affines contrôlent
   la taille fournie à `WRTbl`.
5. Vérifier les invariants des récurrences avec l’ordre incluant le flag et
   préserver le contrat d’exécution numérique.

Les domaines d’entrée restent des hypothèses ou des obligations imposées par
l’hôte. L’analyse des entrées audio et des échantillons soundfile suppose le
domaine fini déclaré ; les widgets supposent des bornes finies et ordonnées,
un initial dans la plage et un pas fini non négatif. Le flag ne borne pas les
zones à l’exécution et ne prouve pas que l’hôte respecte ces déclarations.

Les garanties de la libm cible, les réassociations, FTZ/DAZ, la précision fine du
`lsb`, et les chemins quad/virgule fixe restent les limites documentées du
contrat. Aucun bornage automatique d’un indice dangereux n’est présenté comme
une réparation du programme. Le [corpus de refus](tests/faust/REJECTIONS.md)
reste le témoin de la politique de validation attendue.
