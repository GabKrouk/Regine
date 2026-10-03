# Régine 2022 → Régine 2026 : analyse du code et ce qui a changé

*Analyse du « Regine Pack 2022 » (scripts `Regine install preinstall getReady.R` et `Regine Run.R`, conservés tels quels
dans [`R/legacy_2022/`](R/legacy_2022/)). Tous les chiffres ci-dessous ont été recalculés sur les données du pack et sont
reproduits par [`tests/test_regine.R`](tests/test_regine.R) et [`tests/engine.test.mjs`](tests/engine.test.mjs).*

## Ce que fait Régine 2022

1. Charge le réseau DAP-seq d'O'Malley et al. (2016) : 3 685 526 lignes « TF → gène cible », 387 TF, 32 605 cibles.
2. Construit pour chaque gène la liste des TF qui le fixent (« memoization »), et compte les connexions de chaque TF.
3. `Regine.plot.gb` : pour chaque taille de kernel k = 2…9, tire 100 kernels de k gènes dans la liste et 100 kernels de k
   gènes au hasard, compte les TF communs à tous les gènes du kernel, et fait un boxplot liste vs hasard. On choisit k à l'œil.
4. `Find.me.direct.regulators.KERNEL.STAT.gb` : tire `rep` kernels de taille k dans la liste et `rep` au hasard, compte
   combien de fois chaque TF fixe tout le kernel (*Real*) ou tout le kernel aléatoire (*Noise*), classe par *Real − Noise*.

## Mon avis

**L'idée est bonne.** Le DAP-seq est très permissif : un TF médian fixe 22 % des gènes du génome, un gène médian est fixé
par 84 TF, et deux gènes pris au hasard partagent déjà ~28 TF. Demander quels TF fixent *tous* les gènes d'un kernel, et
comparer à des kernels aléatoires, est une manière simple et intuitive de faire ressortir les régulateurs communs.
Le boxplot liste vs hasard est aussi un très bon diagnostic global (« ma liste partage-t-elle des régulateurs ? »).

**Mais le code a des défauts pratiques et l'approche statistique peut être nettement améliorée.**

### Problèmes de code

| Problème | Conséquence |
|---|---|
| La « memoization » parcourt les 3,68 M lignes pour chacun des 32 606 gènes (`which(a[,2]==gene)`) | ~10¹¹ comparaisons : « takes a Loooong while ». `split()` fait la même chose en 1 s |
| 836 598 lignes du tableau sont des doublons exacts (23 %) | le compteur de connexions (`grep` sur 3,68 M lignes) surestime les cibles de 181 TF, jusqu'à ×1,85 |
| `grep(gene, …)` est une recherche de motif, pas une égalité | fragile ; et la matrice résultat est en texte, ce qui convertit tout le tableau final en caractères |
| Toutes les fonctions lisent des variables globales (`DAP.Seq.Network`, `DAP.Seq.Network.Fast.list`…) | impossible de réutiliser les fonctions sans relancer tout le script dans le bon ordre |
| `install.packages('parallel')` (paquet de base), `library(devtools)` inutile, `mclapply` | erreurs à l'installation ; aucun parallélisme sous Windows |
| Pas de `set.seed` | deux lancements donnent deux tableaux différents (voir plus bas) |
| Un gène absent du réseau donne `NULL` dans `Reduce(intersect, …)` | tout kernel qui le contient compte 0 régulateur, sans avertissement |
| `sample(gene.list, k)` avec k > taille de la liste | erreur |
| `Regine.plot` : `rbind` dans une boucle, nombres stockés en texte, `log(0) = -Inf` dans `qplot` (obsolète) | lent, avertissements |
| Le message final annonce `DAP.seq.Results.txt` | le fichier s'appelle en réalité `DAP.seq.results.<nom>.size….txt` |

### Problèmes de méthode

1. **On tire au hasard quelque chose qui se calcule exactement.** La probabilité qu'un kernel de k gènes de la liste soit
   entièrement fixé par un TF qui fixe k<sub>t</sub> des n gènes de la liste vaut exactement C(k<sub>t</sub>, k) / C(n, k).
   Idem pour le bruit avec C(K<sub>t</sub>, k) / C(N, k). Le tirage de 2022 n'est qu'une estimation bruitée de ces valeurs
   (vérifié : avec 20 000 tirages, la procédure 2022 converge vers les valeurs exactes à < 1,5 % près).
   Avec `rep = 1000` et k = 6 sur la liste nitrate de 13 gènes (1 716 kernels possibles), **le top 10 change de 2 à 3 TF
   d'un lancement à l'autre**, et NLP7 sort du top 10 dans un lancement sur cinq.
2. **Real − Noise favorise les TF « hubs ».** C'est une différence de fréquences, pas un enrichissement. Comme Real ne dépend
   que de k<sub>t</sub>, un TF qui fixe 41 % du génome (DEAR2) passe devant un TF très spécifique. Sur la liste nitrate,
   **NLP7** (régulateur central connu de la réponse au nitrate) fixe 10 des 13 gènes mais seulement 829 gènes du génome
   (enrichissement ×30, FDR 1·10⁻¹¹) : il n'arrive que 8ᵉ au score, derrière DEAR2, MYB67, RVE7L, HB40… (enrichissement ×1,4–2,1,
   FDR > 0,1). La colonne *TF.out.connections* de 2022 servait à corriger cela à l'œil ; aucune p-value n'était donnée.
3. **Le choix du kernel était manuel**, et le critère « là où la différence est la plus grande » ne marche pas tel quel :
   la différence brute des nombres de TF communs est toujours maximale à k = 2 (38 TF d'écart à k = 2, 8 à k = 5 sur la liste
   nitrate) parce que tout le monde partage beaucoup de TF à k = 2.
4. **Le fond aléatoire est fixé** (tous les gènes du réseau). Pour une liste issue d'un RNA-seq ou d'une puce, le bon fond
   est l'ensemble des gènes exprimés ou mesurés.

Le score et la p-value se complètent : sur les listes d'exemple, le tri par P retrouve les régulateurs attendus pour l'ABA
(AREB3, ABI5, GBF6 : les bZIP ABF/AREB) et pour le nitrate (NLP7), alors que le score Régine retrouve MP/ARF5 en tête pour
l'auxine. Les deux classements sont donc proposés.

## Ce que fait Régine 2026

* **Valeurs exactes** de Real, Noise et Real − Noise (en fréquences de kernels, plus de tirage, résultat identique à chaque
  lancement), plus l'**enrichissement** (×), la **p-value hypergéométrique** et la **FDR** (Benjamini–Hochberg sur les 387 TF).
  Classement par FDR par défaut, ou par score Régine.
* **Taille de kernel automatique** : pour k = 2…10, kernels tirés dans la liste et dans le fond (10 000 par taille sur le web,
  graine fixée, donc reproductible). La séparation des deux distributions est mesurée par l'**AUC** : la probabilité qu'un
  kernel de la liste partage plus de TF qu'un kernel aléatoire (0,5 = aucune différence). Le k retenu est celui de l'AUC
  maximale. C'est la traduction chiffrée de « là où les boxplots violet et bleu sont le plus séparés ». Sur les exemples,
  le k automatique tombe entre 4 et 6 (nitrate 2022 : k = 5, AUC 0,88 ; la valeur choisie à l'œil en 2022 était 6).
* **Annotations des TF** : symbole, nom UniProt, famille, description courte Araport11, résumé TAIR et fonction UniProt
  (TAIR / Araport11 / UniProtKB, mêmes sources que GeneCloud 2026), au lieu du seul numéro AGI. Pour chaque TF, la liste
  des gènes de votre liste qu'il fixe (avec leurs symboles).
* **Fond au choix** : gènes du réseau DAP-seq (comportement 2022), gènes codant des protéines, ou votre propre liste.
* **Réseau dédoublonné** (2 848 928 liens distincts) ; gènes non reconnus ou sans pic DAP-seq signalés.
* **Code R réécrit** ([`R/regine.R`](R/regine.R)) : R de base uniquement, sans variables globales, chargement des données et analyse en quelques secondes,
  sorties PDF + TSV. Il donne les mêmes tableaux que la version web.
* **Version web** ([gabkrouk.github.io/Regine](https://gabkrouk.github.io/Regine/)) : tout tourne dans le navigateur,
  la liste n'est envoyée nulle part.

## Limites qui restent

* Le DAP-seq est *in vitro* (ADN génomique nu, sans chromatine ni partenaires) : un TF candidat reste une hypothèse.
* 387 TF seulement (ceux de la table 2022) ; l'attribution pic → gène est celle de la table 2022.
* La famille des TF est déduite des annotations par des règles textuelles (un seul TF sans famille : AT3G10030). Elle est
  écrite dans [`data/tf_annotation.tsv`](data/tf_annotation.tsv) et peut être corrigée à la main.
* L'AUC est estimée par tirage (±0,003 avec 10 000 tirages) : quand le profil est plat, deux tailles voisines peuvent être
  quasi ex æquo. La graine fixe rend le choix reproductible ; le graphique montre l'AUC de chaque taille.
