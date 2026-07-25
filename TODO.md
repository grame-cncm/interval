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

Refaire le meme raisonnement que pour Mod, en prenant en compte les changements dûs à tie to even.
**Note 2026:** Mod a maintenant un chemin ENTIER (sémantique C : signe de x, magnitude
< |y|, identité si |x| < mmin) — voir README et PRECISION.md ; Rem n'a pas encore
d'équivalent.

## Fonctions trigonométriques et hyperboliques

Cos, Sin, Tan, Sinh, Tanh: Refaire le raisonnement des Taylor fallback.

## Neg et sub

Prendre en compte le wrapping sur le meme principe que pour Add (`-INT_MIN > INT_MAX`).
**Toujours ouvert (2026)** : ni `intervalNeg.cpp` ni `intervalSub.cpp` n'ont de
branche entière avec wrap.

## Probablement d'autres fonctions

Faire une passe sur toutes les primitives pour vérifier que tout est en ordre.

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
