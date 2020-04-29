c----------TEMPLATE TO BE READ  AND MODIFIEDBY Tb_rough_vs_freq.m ------
-                        PARAMETRES SATELLITE                          -
c-----------------------------------------------------------------------
-- Longueur d onde du satellite (m) (0.015503876, 0.0081081081)
XXlambda
-- Nombre d onde de coupure (rad/m)
XXKd
9999
-- Angles d incidence du satellite (deg)
XXtht
9999 
c-----------------------------------------------------------------------
-                              VENT                                    -
c-----------------------------------------------------------------------
-- force du Vent (m/s)
XXU10
9999
-- Altitude du vent (m)
XXUalt
-- Modele de Coefficient de trainee(y : Cardone (69), c : Charnock (55), d:Donelan(93))
y
c-----------------------------------------------------------------------
-                              OCEAN                                   -
c-----------------------------------------------------------------------
-- Temperatures de surface de l ocean (deg C)
XXSST
9999
-- Salinité de surface (Psu)
XXSSS
9999
-- Stabilité (Tair - Teau)
XXstab
9999
c-----------------------------------------------------------------------
-                              SORTIES                                 -
c-----------------------------------------------------------------------
-- Fichier de sortie avec écume (aucun => pas de sortie avec écume)
XXFout1
-- Fichier de sortie sans écume (aucun => pas de sortie sans écume)
XXFout2
c-----------------------------------------------------------------------
-                           ETAT DE MER                                -
c-----------------------------------------------------------------------
-- Calcul des Variances (c : Cox & Munk 55, s: Spectre)
s
-- Amplitude du spectre de Durden & Vesecky (0.008 Yueh et 0.004 D & V)
XXSpecAmp
-- Modele de Spectre (e: Elfouhaily, d:Durden & Vesecky)
d
-- Inverse de l''age des vagues 
0.84
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
5.
-- lambdamin
0.
-- lambdamax
0.
-- iVal
0.
-- Type de modele (2 : 2 echelles, p : petite echelles, g : grandes echelles)
2
-- Nombre d integrations sur pente upwind
200
-- Nombre d integrations sur pente xwind
200
-- Limite sup d integration sur pente (X x variance)
5.
