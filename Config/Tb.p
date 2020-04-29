c-----------------------------------------------------------------------
-                        PARAMETRES SATELLITE                          -
c-----------------------------------------------------------------------
-- Longueur d onde du satellite (m) (0.015503876, 0.0081081081)
0.2123
-- Nombre d onde de coupure (rad/m)
4.4842
9999
-- Angles d incidence du satellite (deg)
0.
5.
10.
15.
20.
25.
30.
35.
40.
45.
50.
55.
60.
9999 
c-----------------------------------------------------------------------
-                              VENT                                    -
c-----------------------------------------------------------------------
-- force du Vent (m/s)
4.
5.
6.
9999
-- Altitude du vent (m)
10.
-- Modele de Coefficient de trainee(y : Cardone (69), c : Charnock (55), d:Donelan(93))
y
c-----------------------------------------------------------------------
-                              OCEAN                                   -
c-----------------------------------------------------------------------
-- Temperatures de surface de l ocean (deg C)
15.
9999
-- Salinité de surface (Psu)
36.
9999
-- Stabilité (Tair - Teau)
0.
0.
9999
c-----------------------------------------------------------------------
-                              SORTIES                                 -
c-----------------------------------------------------------------------
-- Fichier de sortie avec écume (aucun => pas de sortie avec écume)
aucun
-- Fichier de sortie sans écume (aucun => pas de sortie sans écume)
test.dat
c-----------------------------------------------------------------------
-                           ETAT DE MER                                -
c-----------------------------------------------------------------------
-- Calcul des Variances (c : Cox & Munk 55, s: Spectre)
s
-- Amplitude du spectre de Durden & Vesecky (0.008 Yueh et 0.004 D & V)
0.008
-- Modele de Spectre (e: Elfouhaily, d:Durden & Vesecky)
d
-- Inverse de l''age des vagues 
1.6
c-----------------------------------------------------------------------
-                              ECUME                                   -
c-----------------------------------------------------------------------
-- Couverture d ecume (Monahan (86) : M1, M2 (moindre carres), Monahan & Lu (90): M3 (actif + passif),M4 (actif), WISE2001 : M5)
M1
-- Modele d emissivite de l ecume (S : Stogryn 72)
s
c-----------------------------------------------------------------------
-                           PERMITTIVITE                               -
c-----------------------------------------------------------------------
-- Modele de constante dielectrique (e: ellison, k:klein & swift)
k
-- Partie imaginaire de constante dielectrique > 0 (O/N)
o
c-----------------------------------------------------------------------
-                             HOULE                                    -
c-----------------------------------------------------------------------
-- Presence de houle (O/N)
N
-- Nombre de point de tabulation de la houle (direction du vent)
100
-- Nombre de point de tabulation de la houle (normale au vent)
100
-- RMS des hauteurs de la houle (m)
1.
-- Largeur à mi puissance de la houle (direction du vent) (m-1)
0.0025
-- Largeur à mi puissance de la houle (normale du vent) (m-1)
0.0025
-- pic densite spectrale houle (direction du vent) (m-1)
0.
-- pic densite spectrale houle (normale du vent) (m-1)
0.0314
c-----------------------------------------------------------------------
-                        TABULATION DIFFUSION                          -
c-----------------------------------------------------------------------
-- Valeur de theta max pour la tabulation de la diffusion (deg), pas
90. 
10.
-- lambdamin
0.
-- lambdamax
0.
-- iVal
0.
-- Type de modele (2 : 2 echelles, p : petite echelles, g : grandes echelles)
p
-- Nombre d integrations sur pente upwind
200
-- Nombre d integrations sur pente xwind
200
-- Limite sup d integration sur pente (X x variance)
5.
