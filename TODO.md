# Tests d'intervalles prenant en compte la précision

Deux possibilités:
* soit on impose que la précision soit exactement la meme
* soit on ne l'impose pas

Ou est utilisé == ?
* dans <=, donc ca a son importance
i<=j <=> reunion(i, j) == j
reunion prend en compte la précision: le résultat a la précision la plus fine des deux opérandes

i<=j pourrait prendre la précision en considération, et considérer que la précision de i doit etre plus grossière que celle de j (toutes les valeurs de i sont présentes dans j)

Ou est utilisé <=?
Dans `check.cpp > analyzeUnaryMethod`, mais en conjonction avec un test sur les précisions (à l'extérieur).

**Conclusion:** Il vaut peut-etre mieux que la précision reste en-dehors des tests, quitte à la rajouter après.
**Note 2026:** `checkExact` compare désormais bornes ET lsb strictement ; c'est la
forme utilisée par la suite de régression (voir TESTING.md).

# Erreurs dans les tests

## Singletons

La précision mesurée d'un singleton est de `INT_MAX`, car il faut au moins deux points dans un intervalle pour mesurer correctement une précision.

**Solution:** Ne pas comparer la précision lorsqu'on est sur un singleton.

# Fonctions d'intervalles à compléter

## Fonctions bitwise

And, Or, Xor: il faut expliquer comment on obtient les bornes des intervalles (pour Yann car c'est lui qui a écrit ce code).

## Inv

Expliquer les expériences faites et pourquoi ca n'était pas concluant.

## Mul

Prouver l'optimalité de la précision $l_x + l_y$ par rapport à la pseudo-injectivité.

## Pow

Expliquer l'implémentation de la précision pour les puissances entières.
A reprendre du reliquat de texte à la fin.

## Rem

**Fait (2026-07):** Remainder (IEEE, quotient ties-to-even) n'est plus un placeholder :
résultat dans [−|y|/2, +|y|/2], magnitude bornée par |x|, identité quand |x| < mmin/2.
Reste ouvert : l'analyse de PRÉCISION fine (le raisonnement pseudo-injectivité de Mod
adapté au tie-to-even) — la borne actuelle est exacte sur les bornes, standard sur le
lsb (min des deux).

## Fonctions trigonométriques et hyperboliques

Cos, Sin, Tan, Sinh, Tanh: Refaire le raisonnement des Taylor fallback.

## Neg et sub

**Fait (2026-07)** : branches entières en miroir de l'idiome d'Add — plein int32 en cas
de chevauchement de frontière, wrap cohérent si entièrement au-delà, exact sinon.

## Passe générale de conformité computationnelle

**Faite (2026-07)** sur les opérations entières : Add/Sub/Neg/Mod/Mul/Pow(iPow)/
IntCast/Min/Max/ARsh conformes ; corrigés : Abs (INT_MIN restait atteignable —
abs(INT_MIN) wrappe), Lsh (le facteur 2^k au lsb flottant court-circuitait
l'élargissement entier de Mul), LRsh (k pouvant être 0 laisse passer l'opérande
négatif ; borne 2^(32−k)−1 au lieu de +inf). La question de plateforme (saturation
ARM vs INT_MIN x86 pour un cast hors-plage) est RÉSOLUE sans choix (2026-07) : le
compilateur émettant un cast nu, IntCast modélise l'UNION des deux comportements —
plage int32 pleine dès qu'un dépassement est possible. Si un jour on veut la
portabilité des VALEURS (même sortie sur toute plateforme), c'est une décision de
codegen (cast clampé, éventuellement derrière une option) — chantier séparé.

# Couche affine

- Documenter le raffinement possible des opérations MIXTES (autres que Mul/Div par un
  opérande de pente nulle).
- Étudier une récupération dépendante de T pour les rampes réinitialisées (le widening
  escalade aujourd'hui vers le top du monde au lieu de proposer une pente rattachée à
  la réinitialisation).

# Documentation

- Méthodes de test avec les différents pitfall/edge cases (singletons, ensembles vides etc)
- Backwards propagation
- Fuzzy precision
