# Breitling 100/24 – Trajectory Optimization

Résolution et optimisation du problème de la Coupe Breitling 100/24 : recherche de la trajectoire aérienne minimale reliant deux aérodromes ($s$ et $t$) en visitant au moins $V_{min}$ aérodromes et en couvrant l'ensemble des régions sous une contrainte de portée maximale $R_{max}$.

Projet réalisé dans le cadre du cours de Recherche Opérationnelle (CORO) à l'ENSIIE par moi ( Ziad Abahamou ) et Yassine Boicra.

## 🛠️ Technologies
* **Langage :** Julia
* **Méthodes :** PLNE / MILP & Heuristiques

## 💡 Approches Implémentées

### 1. Modèles Exacts (PLNE)
* **Formulation MTZ :** Modèle polynomial avec variables d'ordre $O(n^2).
* **Élimination itérative de sous-tours :** Détection et ajout dynamique de coupes (Lazy Constraints).

### 2. Heuristiques & Matheuristiques
* **Glouton + BFS :** Construction progressive avec recherche en largeur pour éviter les blocages locaux.
* **Recherche Locale :** Optimisation par mouvements `2-opt` et `Or-opt`.
* **Matheuristique Branch-and-Cut :** Résolution MIP bornée en temps initialisée par la solution gloutonne.

## 🚀 Exécution

```bash
# Cloner le dépôt
git clone [https://github.com/ABAHAMOU/breitling-trajectory-optimization.git](https://github.com/ABAHAMOU/breitling-trajectory-optimization.git)
cd breitling-trajectory-optimization

# Lancer la résolution avec Julia
julia main.jl
