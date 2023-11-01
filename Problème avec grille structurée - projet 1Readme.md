# Problème avec grille structurée

Notre projet a pour objectif de résoudre numériquement un problème au dérivées partiels en utilisant des méthodes itératives, Jacobi et Gauss-Seidel puis des méthodes parallélisées. Nous débutons par la création d'un code séquentiel dit de Jacobi, suivi de sa parallélisation et nous évaluons les performances des deux algorithmes sur une station de travail de l'ENSTA. Parallèlement, nous développons également un code séquentiel pour la méthode de Gauss-Seidel. L'objectif est de comparer les performances des différentes méthodes (après les avoir validées). Cette phase initiale est cruciale pour établir une base solide pour la résolution ultérieure du problème et l'optimisation des méthodes numériques.

Une fois cela fait nous tester notre code et calculons notre résultat via le mésocentre Cholesky (Ecole Polytehcnique). Cela nous permet de calculer de façon plus pertinente l'efficacité de notre parallélisation.

Finalement, nous construisons également une version parallélisée de l'algorithme de Gauss-Siedel via une méthode dite de "Rouge et Noir".

# Dans ce dossier vous trouverez

- Le code "Jacobi.cpp" continent le code suivant la méthode de Jacobi ;
- Le code "Jacobi_para.cpp" continent le code paralléisé de la méthode de Jacobi ;
- Le code "Gauss-Siedel.cpp" continent le code suivant la méthode de Gauss-Siedel ;
- Le code "Gauss-Siedel_para.cpp" continent le code paralléisé de la méthode de Gauss-Siedel ;
- Notre rapport de projet au format pdf

