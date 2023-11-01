# Problème stationnaire, grille non structurée

L'objectif de ce projet est de compléter et d'améliorer un code parrallèle pour un calcul d'élements fini. Un code C++ de résolution d'éléments fini nous est donné et dans un premier temps nous devions ajouter le calcul des residu dans la méthode de Jacobi ainsi qu'un calcul de l'erreur L2 par rapport à une solution exact connue. Une fois ce code validé nous devions construire un nouveau solver utilisant la méthode de gradient conjugué puis le parraléliser.

# Dans ce dossier vous trouverez les fichiers suivants

- main.cpp
- header.hpp
- parallel.cpp (qui regroupe les fonctions pour la parallélisation du code)
- solver.cpp (le solver via la méthode de Jacobi)
- solver_cgm.cpp (le solver via la du gradient conjugué)
- fonction.cpp (qui regroupe les fonctions que nous avons codées : norm_2, erreur_L2, calcul_resid)
- Nos diapositives résumant notre travail
