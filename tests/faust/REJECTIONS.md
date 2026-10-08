---
title: Programmes Faust à refuser
subtitle: Validation statique des indices de table
date: 8 octobre 2026
document-style: article-a4
language: fr
---

::: toc+
- **Indice d’écriture négatif** — isoler une erreur d’indice sans débordement du compteur.
- **Modulo hors de la boucle et wrapping** — couvrir tout int32 avant de borner le reste signé.
- **Résultat attendu et reproduction** — exiger un diagnostic, et vérifier un programme explicitement corrigé.
- **Portée du constat** — distinguer la bibliothèque d’intervalles et son utilisation par Faust.
:::

# Indice d’écriture négatif

Le témoin [rwtable-negative-write-index.dsp](reject/rwtable-negative-write-index.dsp)
utilise un compteur dont le modulo est placé **dans** la boucle :

```faust
i = ((+(1) : %(200)) ~ _) : -(50);
process = rwtable(100, 0.0, i, _, max(0, min(i, 99)));
```

Le compteur reste dans `[0, 199]` et ne déborde pas. Faust calcule
`i ∈ [-50, 149]`. Dès le premier échantillon, `i = -49`. La table
n’accepte que les indices de `[0, 99]`.

L’indice de lecture est explicitement borné par l’auteur. L’indice d’écriture
ne l’est pas : le compilateur doit refuser ce programme. Le bornage automatique
de l’écriture masquerait l’erreur et modifierait le comportement demandé.

# Modulo hors de la boucle et wrapping

Le second témoin,
[rwtable-wrapping-write-index.dsp](reject/rwtable-wrapping-write-index.dsp),
place le modulo après la boucle :

```faust
i = (+(1) ~ _) : %(200) : -(50);
process = rwtable(100, 0.0, i, _, max(0, min(i, 99)));
```

Cette écriture ne borne pas l’état récursif. Sous le contrat d’arithmétique
entière int32 avec wrapping, le compteur atteint `2147483647`, puis passe à
`-2147483648`. Sur une exécution suffisamment longue, il parcourt tout int32.
Son intervalle conservateur doit donc couvrir
`[-2147483648, 2147483647]`.

Il faut distinguer cet intervalle de celui de l’indice final. Le reste entier
`% 200` a le signe du dividende : après wrapping, il peut être négatif.
La propagation donne :

| Étape | Intervalle |
| :--- | :--- |
| État du compteur `(+(1) ~ _)` | `[-2147483648, 2147483647]` |
| Reste entier `% 200` | `[-199, 199]` |
| Indice après `−50` | `[-249, 149]` |

Ce programme doit aussi être refusé. Il vérifie que l’analyse ne confond pas un
modulo extérieur avec une borne de l’état récursif, et qu’elle conserve les
conséquences du wrapping. L’erreur d’accès existe déjà au premier échantillon,
où l’indice vaut `-49` ; la borne supplémentaire `-249` concerne les valeurs
possibles après le wrapping.

# Résultat attendu et reproduction

[run-rejections.py](run-rejections.py) vérifie le refus des deux témoins avec les backends `cpp`
et `ocpp`, en simple et double précision. Il vérifie aussi l’acceptation du
[témoin valide](accept/rwtable-explicitly-bounded-indices.dsp), où l’auteur borne
explicitement les deux indices. Une interruption ou un crash du compilateur
ne vaut pas un refus avec diagnostic.

```sh
python3 tests/faust/run-rejections.py --faust /chemin/vers/faust
```

Avec le binaire d’intégration `669e8de5b-modified`, le témoin valide est accepté
dans les quatre configurations. Les deux témoins dangereux sont également acceptés dans
les quatre configurations : **les huit vérifications de refus échouent actuellement** et le
pilote retourne un statut non nul. Il n’est pas transformé en succès attendu.
Ce test dépend d’un compilateur Faust externe et reste séparé du CTest de la
bibliothèque autonome.

# Portée du constat

L’analyse de l’indice contient bien les valeurs dangereuses. Le même indice,
utilisé dans `process = _ @ i`, provoque un refus de compilation signalant
`interval(-50,149,0)`. Ce témoin ne démontre donc pas une erreur de bornes de la
bibliothèque : il démontre une validation insuffisante des indices de table par
Faust. La validation recherchée exige un intervalle entièrement contenu dans
`[0, N − 1]` ; une impossibilité de certifier cette inclusion doit rester visible.

Pour la bibliothèque, les tests achevés passent en float et double avec Clang et
GCC, avec les oracles MPFR facultatifs, ainsi que sous sanitizers et en
WebAssembly. Aucun contre-exemple restant n’a été identifié dans ces tests.
Cela ne constitue pas une certification complète. Les limites déjà documentées
dans le [README](../../README.md) restent à traiter :

- NaN n’est pas représenté par un indicateur de valeur possible séparé. Un
  intervalle numérique ne suffit donc pas à certifier qu’un résultat est défini.
- `IntCast` renvoie tout int32 lorsque son entrée peut sortir de cette plage.
  Cette enveloppe couvre les comportements matériels considérés, mais ne rend
  pas défini un cast C++ hors plage. La validité de la conversion doit être
  établie séparément avant de certifier un indice ; il manque une représentation
  explicite de cette impossibilité de certifier.
- La compensation libm suppose au plus deux pas représentables depuis le
  résultat correctement arrondi ; aucune garantie universelle par cible n’est
  fournie.
- Les invariants des récurrences doivent être vérifiés au niveau du solveur.
  Les règles locales d’arrondi ne suffisent pas à établir ces invariants.
- Les chemins quad et virgule fixe et la précision fine du `lsb` n’ont pas reçu
  la même certification que les bornes float et double.

Le contrat d’exécution, les domaines d’entrée et les transformations effectuées
par le compilateur doivent aussi correspondre à ce que l’analyse modélise.
