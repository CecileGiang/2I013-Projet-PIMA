# Projet-PIMA_2I013
(Projet PIMA 2017-2018) Evaluer vite et bien: de la robotique aux polynômes

**Encadrant:** Mohab Safey El Din

**Etudiants:** GIANG Cécile, OSUNA VARGAS Ana Pamela,YE Weijie

**Titre: Evaluer vite et bien: de la robotique aux polynômes**

**_Description:_** Ce projet est motivé par une application et robotique qui faitl'objet d'une collaboration entre l'UPMC (Sorbonne Université) et l'université J.Kepler en Autriche: la planification de trajectoire de robots dont les propriétés mécaniques sont particulières. Après un long processus de modélisation, on est ramené à étudierle nombre de composantes connexes de courbes planes (définies par l'annulation de polynômes en deux variables).
Cette étape donne lieu à des calculs intensifs et très coûteux. Ils nécessitent  d'évaluer des polynômes en une variable de degré gigantesque (autour de 100 000) et dont les tailles de coefficients sont tout aussi imposantes (de l'ordre du million de bits).
Il s'agit donc de développer des techniques d'algorithmique fine combinées à du calcul haute-performance, et des techniques de programmation astucieuses pour pouvoir traiter des données de cette taille. Concrètement, nous développerons une algorithmique _diviser-pour-régner_ asymptotiquement optimale combinée à des astuces permettant de limiter la précisieon des calculs (néanmoins effectuées en arithmétique multi-précision). Enfin, l'algorithme sera parallélisé; on visera une utilisation optimale des threads disponibles.

**Langage utilisé:** C

**Outils:** Bibliothèque GNU GMP pour la manipulation de grands entiers
