# SimulationPhotocatalysis

Ce dossier contient :
- le fichier de données obtenu en opérant le photoréacteur en batch : _**"Test.xlsx"**_
- le fichier de données obtenu en opérant le photoréacteur en continu : _**"TestContinu.xlsx"**_
- les fichiers d'irradiation pour les villes de Bordeaux, Nantes et Montpellier
- les différentes fonctions matlab permettant de :
        - résoudre les bilans de matière sur le réacteur avec les différentes lois cinétiques :_**"ResolBilan.m"**_ en batch et _**"ResolBilan_SimulationContinu.m"**_ en continu
        - optimiser le calcul des paramètre de chaque loi cinétique : _**"Optimisation.m"**_
        - calculer la durée du plateau initiale à partir des données expérimentales en batch : _**"InitialShoulder.m"**_
        - moduler le débit d'alimentation pour le fonctionnement du réacteur en continu : _**"ResolBilan_SimulationDebit.m"**_
        - calculer les vitesse d'inactivation moyenne dans le réacteur : _**"VitesseInactivation_Moyenne.m"**_
- le MAIN contenant le corps du code et permettant de faire tourner les simulations à partir des différentes donctions et données expérimetales : _**"CodeSimulationArticle_I_lisse.mlx"**_
